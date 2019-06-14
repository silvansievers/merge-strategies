#include "merge_and_shrink_algorithm.h"

#include "distances.h"
#include "factored_transition_system.h"
#include "fts_factory.h"
#include "label_reduction.h"
#include "labels.h"
#include "merge_and_shrink_representation.h"
#include "merge_strategy.h"
#include "merge_strategy_factory.h"
#include "shrink_strategy.h"
#include "transition_system.h"
#include "types.h"
#include "utils.h"

#include "../options/option_parser.h"
#include "../options/options.h"

#include "../task_utils/task_properties.h"

#include "../utils/countdown_timer.h"
#include "../utils/logging.h"
#include "../utils/markup.h"
#include "../utils/math.h"
#include "../utils/system.h"
#include "../utils/timer.h"

#include <cassert>
#include <iostream>
#include <set>
#include <string>
#include <utility>
#include <vector>

using namespace std;
using options::Bounds;
using options::OptionParser;
using options::Options;
using utils::ExitCode;

namespace merge_and_shrink {
static void log_progress(const utils::Timer &timer, string msg) {
    cout << "M&S algorithm timer: " << timer << " (" << msg << ")" << endl;
}

MergeAndShrinkAlgorithm::MergeAndShrinkAlgorithm(const Options &opts) :
    merge_strategy_factory(opts.get<shared_ptr<MergeStrategyFactory>>("merge_strategy")),
    shrink_strategy(opts.get<shared_ptr<ShrinkStrategy>>("shrink_strategy")),
    label_reduction(opts.get<shared_ptr<LabelReduction>>("label_reduction", nullptr)),
    max_states(opts.get<int>("max_states")),
    max_states_before_merge(opts.get<int>("max_states_before_merge")),
    shrink_threshold_before_merge(opts.get<int>("threshold_before_merge")),
    prune_unreachable_states(opts.get<bool>("prune_unreachable_states")),
    prune_irrelevant_states(opts.get<bool>("prune_irrelevant_states")),
    pruning_as_abstraction(opts.get<bool>("pruning_as_abstraction")),
    verbosity(static_cast<utils::Verbosity>(opts.get_enum("verbosity"))),
    main_loop_max_time(opts.get<double>("main_loop_max_time")),
    num_transitions_to_abort(opts.get<int>("num_transitions_to_abort")),
    num_transitions_to_exclude(opts.get<int>("num_transitions_to_exclude")),
    starting_peak_memory(0) {
    assert(max_states_before_merge > 0);
    assert(max_states >= max_states_before_merge);
    assert(shrink_threshold_before_merge <= max_states_before_merge);
}

void MergeAndShrinkAlgorithm::report_peak_memory_delta(bool final) const {
    if (final)
        cout << "Final";
    else
        cout << "Current";
    cout << " peak memory increase of merge-and-shrink algorithm: "
         << utils::get_peak_memory_in_kb() - starting_peak_memory << " KB"
         << endl;
}

void MergeAndShrinkAlgorithm::dump_options() const {
    if (verbosity >= utils::Verbosity::VERBOSE) {
        if (merge_strategy_factory) { // deleted after merge strategy extraction
            merge_strategy_factory->dump_options();
            cout << endl;
        }

        cout << "Options related to size limits and shrinking: " << endl;
        cout << "Transition system size limit: " << max_states << endl
             << "Transition system size limit right before merge: "
             << max_states_before_merge << endl;
        cout << "Threshold to trigger shrinking right before merge: "
             << shrink_threshold_before_merge << endl;
        cout << endl;

        shrink_strategy->dump_options();
        cout << endl;

        if (label_reduction) {
            label_reduction->dump_options();
        } else {
            cout << "Label reduction disabled" << endl;
        }
        cout << endl;
    }
}

void MergeAndShrinkAlgorithm::warn_on_unusual_options() const {
    string dashes(79, '=');
    if (!label_reduction) {
        cout << dashes << endl
             << "WARNING! You did not enable label reduction.\nThis may "
            "drastically reduce the performance of merge-and-shrink!"
             << endl << dashes << endl;
    } else if (label_reduction->reduce_before_merging() && label_reduction->reduce_before_shrinking()) {
        cout << dashes << endl
             << "WARNING! You set label reduction to be applied twice in each merge-and-shrink\n"
            "iteration, both before shrinking and merging. This double computation effort\n"
            "does not pay off for most configurations!"
             << endl << dashes << endl;
    } else {
        if (label_reduction->reduce_before_shrinking() &&
            (shrink_strategy->get_name() == "f-preserving"
             || shrink_strategy->get_name() == "random")) {
            cout << dashes << endl
                 << "WARNING! Bucket-based shrink strategies such as f-preserving random perform\n"
                "best if used with label reduction before merging, not before shrinking!"
                 << endl << dashes << endl;
        }
        if (label_reduction->reduce_before_merging() &&
            shrink_strategy->get_name() == "bisimulation") {
            cout << dashes << endl
                 << "WARNING! Shrinking based on bisimulation performs best if used with label\n"
                "reduction before shrinking, not before merging!"
                 << endl << dashes << endl;
        }
    }

    if (!prune_unreachable_states || !prune_irrelevant_states) {
        cout << dashes << endl
             << "WARNING! Pruning is (partially) turned off!\nThis may "
            "drastically reduce the performance of merge-and-shrink!"
             << endl << dashes << endl;
    }
}

bool MergeAndShrinkAlgorithm::ran_out_of_time(
    const utils::CountdownTimer &timer) const {
    if (timer.is_expired()) {
        if (verbosity >= utils::Verbosity::NORMAL) {
            cout << "Ran out of time, stopping computation." << endl;
            cout << endl;
        }
        return true;
    }
    return false;
}

bool MergeAndShrinkAlgorithm::too_many_transitions(const FactoredTransitionSystem &fts, int index) const {
    int num_transitions = fts.get_transition_system(index).compute_total_transitions();
    if (num_transitions > num_transitions_to_abort) {
        if (verbosity >= utils::Verbosity::NORMAL) {
            cout << "Factor has too many transitions, stopping computation."
                 << endl;
            cout << endl;
        }
        return true;
    }
    return false;
}

bool MergeAndShrinkAlgorithm::too_many_transitions(const FactoredTransitionSystem &fts) const {
    for (int index = 0; index < fts.get_size(); ++index) {
        if (fts.is_active(index)) {
            if (too_many_transitions(fts, index)) {
                return true;
            }
        }
    }
    return false;
}

bool MergeAndShrinkAlgorithm::exclude_if_too_many_transitions() const {
    return num_transitions_to_exclude != INF;
}

void MergeAndShrinkAlgorithm::main_loop(
    FactoredTransitionSystem &fts,
    const TaskProxy &task_proxy) {
    utils::CountdownTimer timer(main_loop_max_time);
    if (verbosity >= utils::Verbosity::NORMAL) {
        cout << "Starting main loop ";
        if (main_loop_max_time == numeric_limits<double>::infinity()) {
            cout << "without a time limit." << endl;
        } else {
            cout << "with a time limit of "
                 << main_loop_max_time << "s." << endl;
        }
    }
    int maximum_intermediate_size = 0;
    int maximum_transitions_size = 0;
    for (int i = 0; i < fts.get_size(); ++i) {
        int size = fts.get_transition_system(i).get_size();
        if (size > maximum_intermediate_size) {
            maximum_intermediate_size = size;
        }
        int num_trans = fts.get_transition_system(i).compute_total_transitions();
        if (num_trans > maximum_transitions_size) {
            maximum_transitions_size = num_trans;
        }
    }

    pair<int, int> score_based_merging_tiebreaking;
    vector<int> init_hvalue_increase;
    vector<int> remaining_labels;
    remaining_labels.push_back(fts.get_labels().compute_number_active_labels());
    bool still_perfect = true;
    vector<pair<int, int>> merge_order;
    vector<double> relative_pruning_per_iteration;
    int num_attempts_merging_for_symmetries = 0;
    int num_imperfect_shrinking_merging_for_symmetries = 0;
    int num_pruning_merging_for_symmetries = 0;
    int num_failed_merging_for_symmetries = 0;
    bool merging_for_symmetries = true;
    bool currently_shrink_perfect_for_symmetries = true;
    bool currently_prune_perfect_for_symmetries = true;

    if (label_reduction) {
        label_reduction->initialize(task_proxy);
    }
    unique_ptr<MergeStrategy> merge_strategy =
        merge_strategy_factory->compute_merge_strategy(task_proxy, fts);
    merge_strategy_factory = nullptr;

    auto log_main_loop_progress = [&timer](const string &msg) {
            cout << "M&S algorithm main loop timer: "
                 << timer.get_elapsed_time()
                 << " (" << msg << ")" << endl;
        };
    int iteration_counter = 0;
    set<int> allowed_indices;
    while (fts.get_num_active_entries() > 1) {
        // Choose next transition systems to merge
        vector<int> vec_allowed_indices;;
        if (exclude_if_too_many_transitions()) {
            vec_allowed_indices = vector<int>(
                allowed_indices.begin(), allowed_indices.end());
        }
        pair<int, int> merge_indices = merge_strategy->get_next(vec_allowed_indices);
        if (ran_out_of_time(timer)) {
            break;
        }
        if (merge_strategy->ended_merging_for_symmetries()) {
            merging_for_symmetries = false;
            if (!currently_shrink_perfect_for_symmetries) {
                ++num_imperfect_shrinking_merging_for_symmetries;
            }
            if (!currently_prune_perfect_for_symmetries) {
                ++num_pruning_merging_for_symmetries;
            }
            if (!currently_shrink_perfect_for_symmetries ||
                !currently_prune_perfect_for_symmetries) {
                ++num_failed_merging_for_symmetries;
            }
        }
        if (merge_strategy->started_merging_for_symmetries()) {
            ++num_attempts_merging_for_symmetries;
            merging_for_symmetries = true;
            currently_shrink_perfect_for_symmetries = true;
            currently_prune_perfect_for_symmetries = true;
        }
        int merge_index1 = merge_indices.first;
        int merge_index2 = merge_indices.second;
        assert(merge_index1 != merge_index2);
        merge_order.push_back(merge_indices);
        if (verbosity >= utils::Verbosity::NORMAL) {
            cout << "Next pair of indices: ("
                 << merge_index1 << ", " << merge_index2 << ")" << endl;
            if (verbosity >= utils::Verbosity::VERBOSE) {
                fts.statistics(merge_index1);
                fts.statistics(merge_index2);
            }
            log_main_loop_progress("after computation of next merge");
        }

        // Label reduction (before shrinking)
        if (label_reduction && label_reduction->reduce_before_shrinking()) {
            bool reduced = label_reduction->reduce(merge_indices, fts, verbosity);
            if (verbosity >= utils::Verbosity::NORMAL && reduced) {
                log_main_loop_progress("after label reduction");
            }
            remaining_labels.push_back(fts.get_labels().compute_number_active_labels());
        }

        if (ran_out_of_time(timer)) {
            break;
        }

        // Shrinking
        bool shrunk = shrink_before_merge_step(
            fts,
            merge_index1,
            merge_index2,
            max_states,
            max_states_before_merge,
            shrink_threshold_before_merge,
            *shrink_strategy,
            verbosity);
        if (verbosity >= utils::Verbosity::NORMAL && shrunk) {
            log_main_loop_progress("after shrinking");
        }
        if (merging_for_symmetries && currently_shrink_perfect_for_symmetries && shrunk) {
            currently_shrink_perfect_for_symmetries = false;
        }

        const vector<double> &miss_qualified_states_ratios =
            shrink_strategy->get_miss_qualified_states_ratios();
        int size = miss_qualified_states_ratios.size();
        if (size >= 2 && still_perfect &&
            (miss_qualified_states_ratios[size - 1]
             || miss_qualified_states_ratios[size - 2])) {
            // The test for size >= 2 is to ensure we actually record
            // this kind of statistics -- currently only with bisimulation
            // shrinking.
            cout << "not perfect anymore in iteration " << iteration_counter << endl;
            still_perfect = false;
        }

        if (ran_out_of_time(timer)) {
            break;
        }

        // Label reduction (before merging)
        if (label_reduction && label_reduction->reduce_before_merging()) {
            bool reduced = label_reduction->reduce(merge_indices, fts, verbosity);
            if (verbosity >= utils::Verbosity::NORMAL && reduced) {
                log_main_loop_progress("after label reduction");
            }
            remaining_labels.push_back(fts.get_labels().compute_number_active_labels());
        }

        if (ran_out_of_time(timer)) {
            break;
        }

        int init_dist1 = fts.get_init_state_goal_distance(merge_index1);
        int init_dist2 = fts.get_init_state_goal_distance(merge_index2);

        // Merging
        int merged_index = fts.merge(merge_index1, merge_index2, verbosity);
        int abs_size = fts.get_transition_system(merged_index).get_size();
        if (abs_size > maximum_intermediate_size) {
            maximum_intermediate_size = abs_size;
        }
        int num_trans = fts.get_transition_system(merged_index).compute_total_transitions();
        if (num_trans > maximum_transitions_size) {
            maximum_transitions_size = num_trans;
        }
        int new_init_dist = fts.get_init_state_goal_distance(merged_index);
        int difference = new_init_dist - max(init_dist1, init_dist2);
        init_hvalue_increase.push_back(difference);

        if (verbosity >= utils::Verbosity::NORMAL) {
            if (verbosity >= utils::Verbosity::VERBOSE) {
                fts.statistics(merged_index);
                cout << "Difference of init h values: " << difference << endl;
            }
            log_main_loop_progress("after merging");
        }

        // We do not check for num transitions here but only after pruning
        // to allow recovering a too large product.
        if (ran_out_of_time(timer)) {
            break;
        }

        // Pruning
        if (prune_unreachable_states || prune_irrelevant_states) {
            int old_size = fts.get_transition_system(merged_index).get_size();
            pair<bool, bool> pruned_and_pruned_unreachable = prune_step(
                fts,
                merged_index,
                prune_unreachable_states,
                prune_irrelevant_states,
                pruning_as_abstraction,
                verbosity);
            double new_size = fts.get_transition_system(merged_index).get_size();
            assert(new_size <= old_size);
            relative_pruning_per_iteration.push_back(1 - new_size / static_cast<double>(old_size));
            if (verbosity >= utils::Verbosity::NORMAL && pruned_and_pruned_unreachable.first) {
                if (verbosity >= utils::Verbosity::VERBOSE) {
                    fts.statistics(merged_index);
                }
                log_main_loop_progress("after pruning");
            }
            if (pruned_and_pruned_unreachable.first &&
                pruned_and_pruned_unreachable.second &&
                merging_for_symmetries &&
                currently_prune_perfect_for_symmetries) {
                currently_prune_perfect_for_symmetries = false;
            }
        }

        /*
          NOTE: both the shrink strategy classes and the construction
          of the composite transition system require the input
          transition systems to be non-empty, i.e. the initial state
          not to be pruned/not to be evaluated as infinity.
        */
        if (!fts.is_factor_solvable(merged_index)) {
            if (verbosity >= utils::Verbosity::NORMAL) {
                cout << "Abstract problem is unsolvable, stopping "
                    "computation. " << endl << endl;
            }
            break;
        }

        if (exclude_if_too_many_transitions()) {
            allowed_indices.erase(merge_index1);
            allowed_indices.erase(merge_index2);
            if (num_trans <= num_transitions_to_exclude) {
                allowed_indices.insert(merged_index);
            } else {
                if (verbosity >= utils::Verbosity::NORMAL) {
                    cout << fts.get_transition_system(merged_index).tag()
                         << "too many number of transitions, excluding "
                            "from further consideration." << endl;
                }
            }
            if (allowed_indices.size() <= 1) {
                if (verbosity >= utils::Verbosity::NORMAL) {
                    cout << "Not enough factors remaining with a low enough "
                            "number of transitions, stopping computation."
                         << endl;
                    cout << endl;
                }
                break;
            }
        }

        if (ran_out_of_time(timer) || too_many_transitions(fts, merged_index)) {
            break;
        }

        // End-of-iteration output.
        if (verbosity >= utils::Verbosity::VERBOSE) {
            report_peak_memory_delta();
        }
        if (verbosity >= utils::Verbosity::NORMAL) {
            cout << endl;
        }

        ++iteration_counter;
    }

    cout << "End of merge-and-shrink algorithm, statistics:" << endl;
    cout << "Main loop runtime: " << timer.get_elapsed_time() << endl;
    cout << "Maximum intermediate abstraction size: "
         << maximum_intermediate_size << endl;
    score_based_merging_tiebreaking =
        merge_strategy->get_tiebreaking_statistics();
    cout << "Iterations with merge tiebreaking: "
         << score_based_merging_tiebreaking.first << endl;
    cout << "Total tiebreaking merge candidates: "
         << score_based_merging_tiebreaking.second << endl;
    cout << "Maximum intermediate number of transitions: "
         << maximum_transitions_size << endl;
    cout << "Init h value improvements: " << init_hvalue_increase << endl;
    cout << "Course of label reduction: " << remaining_labels << endl;
    const vector<double> &miss_qualified_states_ratios =
        shrink_strategy->get_miss_qualified_states_ratios();
    cout << "Course of miss qualified states shrinking: "
         << miss_qualified_states_ratios << endl;
    double summed_values = 0;
    for (double value : miss_qualified_states_ratios) {
        summed_values += value;
    }
    size_t number_of_shrinks = miss_qualified_states_ratios.size();
    double average_imperfect_shrinking = 0;
    if (number_of_shrinks) {
        average_imperfect_shrinking = summed_values / static_cast<double>(number_of_shrinks);
    }
    cout << "Average imperfect shrinking: " << average_imperfect_shrinking << endl;
    cout << "Merge order: [";
    bool linear_order = true;
    int next_index = task_proxy.get_variables().size();
    for (size_t i = 0; i < merge_order.size(); ++i) {
        pair<int, int> merge = merge_order[i];
        cout << "(" << merge.first << ", " << merge.second << ")";
        if (i != merge_order.size() - 1) {
            cout << ", ";
        }
        if (linear_order && i != 0) {
            if (merge.first != next_index && merge.second != next_index) {
                linear_order = false;
            }
            ++next_index;
        }
    }
    cout << "]" << endl;
    if (linear_order) {
        cout << "Linear merge order" << endl;
    } else {
         cout << "Non-linear merge order" << endl;
    }
    cout << "Relative pruning per iteration: " << relative_pruning_per_iteration << endl;
    double summed_pruning = 0;
    for (double pruning : relative_pruning_per_iteration) {
        summed_pruning += pruning;
    }
    // If relative_pruning_per_iteration are empty, then because the instance is unsolvable.
    // In this case, we return 0, which is the worst value possible for pruning.
    double average_pruning = 0;
    if (!relative_pruning_per_iteration.empty()) {
        average_pruning =  summed_pruning / static_cast<double>(relative_pruning_per_iteration.size());
    }
    cout << "Average relative pruning: " << average_pruning << endl;

    cout << "Number of attempts to merge for symmetries: "
         << num_attempts_merging_for_symmetries << endl;
    cout << "Number of times non-perfect shrinking interfered merging for symmetries: "
         << num_imperfect_shrinking_merging_for_symmetries << endl;
    cout << "Number of times pruning interfered merging for symmetries: "
         << num_pruning_merging_for_symmetries << endl;
    cout << "Number of times merging for symmetries failed for any reason: "
         << num_failed_merging_for_symmetries << endl;
    cout << endl;

    shrink_strategy = nullptr;
    label_reduction = nullptr;
}

FactoredTransitionSystem MergeAndShrinkAlgorithm::build_factored_transition_system(
    const TaskProxy &task_proxy) {
    if (starting_peak_memory) {
        cerr << "Calling build_factored_transition_system twice is not "
             << "supported!" << endl;
        utils::exit_with(utils::ExitCode::SEARCH_CRITICAL_ERROR);
    }
    starting_peak_memory = utils::get_peak_memory_in_kb();

    utils::Timer timer;
    cout << "Running merge-and-shrink algorithm..." << endl;
    task_properties::verify_no_axioms(task_proxy);
    dump_options();
    warn_on_unusual_options();
    cout << endl;

    const bool compute_init_distances =
        shrink_strategy->requires_init_distances() ||
        merge_strategy_factory->requires_init_distances() ||
        prune_unreachable_states;
    const bool compute_goal_distances =
        shrink_strategy->requires_goal_distances() ||
        merge_strategy_factory->requires_goal_distances() ||
        prune_irrelevant_states;
    FactoredTransitionSystem fts =
        create_factored_transition_system(
            task_proxy,
            compute_init_distances,
            compute_goal_distances,
            verbosity);
    if (verbosity >= utils::Verbosity::NORMAL) {
        log_progress(timer, "after computation of atomic factors");
    }

    /*
      Prune all atomic factors according to the chosen options. Stop early if
      one factor is unsolvable.

      TODO: think about if we can prune already while creating the atomic FTS.
    */
    bool pruned = false;
    bool unsolvable = false;
    for (int index = 0; index < fts.get_size(); ++index) {
        assert(fts.is_active(index));
        if (prune_unreachable_states || prune_irrelevant_states) {
            pair<bool, bool> pruned_and_pruned_unreachable = prune_step(
                fts,
                index,
                prune_unreachable_states,
                prune_irrelevant_states,
                pruning_as_abstraction,
                verbosity);
            pruned = pruned || pruned_and_pruned_unreachable.first;
        }
        if (!fts.is_factor_solvable(index)) {
            cout << "Atomic FTS is unsolvable, stopping computation." << endl;
            unsolvable = true;
            break;
        }
    }
    if (verbosity >= utils::Verbosity::NORMAL) {
        if (pruned) {
            log_progress(timer, "after pruning atomic factors");
        }
        cout << endl;
    }

    if (!unsolvable && !too_many_transitions(fts) && main_loop_max_time > 0) {
        main_loop(fts, task_proxy);
    }
    const bool final = true;
    report_peak_memory_delta(final);
    cout << "Merge-and-shrink algorithm runtime: " << timer << endl;
    cout << endl;
    return fts;
}

void add_merge_and_shrink_algorithm_options_to_parser(OptionParser &parser) {
    // Merge strategy option.
    parser.add_option<shared_ptr<MergeStrategyFactory>>(
        "merge_strategy",
        "See detailed documentation for merge strategies. "
        "We currently recommend SCC-DFP, which can be achieved using "
        "{{{merge_strategy=merge_sccs(order_of_sccs=topological,merge_selector="
        "score_based_filtering(scoring_functions=[goal_relevance,dfp,total_order"
        "]))}}}");

    // Shrink strategy option.
    parser.add_option<shared_ptr<ShrinkStrategy>>(
        "shrink_strategy",
        "See detailed documentation for shrink strategies. "
        "We currently recommend non-greedy shrink_bisimulation, which can be "
        "achieved using {{{shrink_strategy=shrink_bisimulation(greedy=false)}}}");

    // Label reduction option.
    parser.add_option<shared_ptr<LabelReduction>>(
        "label_reduction",
        "See detailed documentation for labels. There is currently only "
        "one 'option' to use label_reduction, which is {{{label_reduction=exact}}} "
        "Also note the interaction with shrink strategies.",
        OptionParser::NONE);

    // Pruning options.
    parser.add_option<bool>(
        "prune_unreachable_states",
        "If true, prune abstract states unreachable from the initial state.",
        "true");
    parser.add_option<bool>(
        "prune_irrelevant_states",
        "If true, prune abstract states from which no goal state can be "
        "reached.",
        "true");
    parser.add_option<bool>(
        "pruning_as_abstraction",
        "If true, perform pruning by not removing pruned states, but mapping "
        "them to two single states representing unreachable and irrelevant "
        "states respectively.",
        "false");

    add_transition_system_size_limit_options_to_parser(parser);

    /*
      silent: no output during construction, only starting and final statistics
      normal: basic output during construction, starting and final statistics
      verbose: full output during construction, starting and final statistics
    */
    utils::add_verbosity_option_to_parser(parser);

    parser.add_option<double>(
        "main_loop_max_time",
        "A limit in seconds on the runtime of the main loop of the algorithm. "
        "If the limit is exceeded, the algorithm terminates, potentially "
        "returning a factored transition system with several factors. Also "
        "note that the time limit is only checked between transformations "
        "of the main loop, but not during, so it can be exceeded if a "
        "transformation is runtime-intense.",
        "infinity",
        Bounds("0.0", "infinity"));
    parser.add_option<int>(
        "num_transitions_to_abort",
        "A limit on the number of transitions of any factor during the "
        "computation. Once this limit is reached, the algorithm terminates, "
        "leaving the chosen partial_mas_method to compute a heuristic from the "
        "set of remaining factors.",
        "infinity",
        Bounds("0", "infinity"));
    parser.add_option<int>(
        "num_transitions_to_exclude",
        "A limit on the number of transitions of any factor during the "
        "computation. Once a factor reaches this limit, it is excluded from "
        "further considerations of the algorithm.",
        "infinity",
        Bounds("0", "infinity"));
}

void add_transition_system_size_limit_options_to_parser(OptionParser &parser) {
    parser.add_option<int>(
        "max_states",
        "maximum transition system size allowed at any time point.",
        "-1",
        Bounds("-1", "infinity"));
    parser.add_option<int>(
        "max_states_before_merge",
        "maximum transition system size allowed for two transition systems "
        "before being merged to form the synchronized product.",
        "-1",
        Bounds("-1", "infinity"));
    parser.add_option<int>(
        "threshold_before_merge",
        "If a transition system, before being merged, surpasses this soft "
        "transition system size limit, the shrink strategy is called to "
        "possibly shrink the transition system.",
        "-1",
        Bounds("-1", "infinity"));
}

void handle_shrink_limit_options_defaults(Options &opts) {
    int max_states = opts.get<int>("max_states");
    int max_states_before_merge = opts.get<int>("max_states_before_merge");
    int threshold = opts.get<int>("threshold_before_merge");

    // If none of the two state limits has been set: set default limit.
    if (max_states == -1 && max_states_before_merge == -1) {
        max_states = 50000;
    }

    // If exactly one of the max_states options has been set, set the other
    // so that it imposes no further limits.
    if (max_states_before_merge == -1) {
        max_states_before_merge = max_states;
    } else if (max_states == -1) {
        int n = max_states_before_merge;
        if (utils::is_product_within_limit(n, n, INF)) {
            max_states = n * n;
        } else {
            max_states = INF;
        }
    }

    if (max_states_before_merge > max_states) {
        cout << "warning: max_states_before_merge exceeds max_states, "
             << "correcting." << endl;
        max_states_before_merge = max_states;
    }

    if (max_states < 1) {
        cerr << "error: transition system size must be at least 1" << endl;
        utils::exit_with(ExitCode::SEARCH_INPUT_ERROR);
    }

    if (max_states_before_merge < 1) {
        cerr << "error: transition system size before merge must be at least 1"
             << endl;
        utils::exit_with(ExitCode::SEARCH_INPUT_ERROR);
    }

    if (threshold == -1) {
        threshold = max_states;
    }
    if (threshold < 1) {
        cerr << "error: threshold must be at least 1" << endl;
        utils::exit_with(ExitCode::SEARCH_INPUT_ERROR);
    }
    if (threshold > max_states) {
        cout << "warning: threshold exceeds max_states, correcting" << endl;
        threshold = max_states;
    }

    opts.set<int>("max_states", max_states);
    opts.set<int>("max_states_before_merge", max_states_before_merge);
    opts.set<int>("threshold_before_merge", threshold);
}
}
