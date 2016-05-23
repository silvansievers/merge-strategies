#include "merge_and_shrink_heuristic.h"

#include "factored_transition_system.h"
#include "fts_factory.h"
#include "label_reduction.h"
#include "labels.h"
#include "merge_strategy.h"
#include "merge_strategy_factory.h"
#include "shrink_strategy.h"
#include "transition_system.h"

#include "../option_parser.h"
#include "../plugin.h"
#include "../task_tools.h"

#include "../utils/logging.h"
#include "../utils/markup.h"
#include "../utils/math.h"
#include "../utils/memory.h"
#include "../utils/system.h"
#include "../utils/timer.h"

#include <cassert>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

using namespace std;
using utils::ExitCode;

namespace merge_and_shrink {
void print_time(const utils::Timer &timer, string text) {
    cout << "t=" << timer << " (" << text << ")" << endl;
}

MergeAndShrinkHeuristic::MergeAndShrinkHeuristic(const Options &opts)
    : Heuristic(opts),
      merge_strategy_factory(opts.get<shared_ptr<MergeStrategyFactory>>("merge_strategy")),
      shrink_strategy(opts.get<shared_ptr<ShrinkStrategy>>("shrink_strategy")),
      label_reduction(nullptr),
      starting_peak_memory(-1),
      max_states(opts.get<int>("max_states")),
      max_states_before_merge(opts.get<int>("max_states_before_merge")),
      shrink_threshold_before_merge(opts.get<int>("threshold_before_merge")),
      debug_transition_systems(opts.get<bool>("debug_transition_systems")),
      fts(nullptr) {
    assert(max_states_before_merge > 0);
    assert(max_states >= max_states_before_merge);
    assert(shrink_threshold_before_merge <= max_states_before_merge);

    if (opts.contains("label_reduction")) {
        label_reduction = opts.get<shared_ptr<LabelReduction>>("label_reduction");
        label_reduction->initialize(task_proxy);
    }

    utils::Timer timer;
    cout << "Initializing merge-and-shrink heuristic..." << endl;
    starting_peak_memory = utils::get_peak_memory_in_kb();
    verify_no_axioms(task_proxy);
    dump_options();
    warn_on_unusual_options();
    cout << endl;

    build_transition_system(timer);
    report_peak_memory_delta(true);
    cout << "Done initializing merge-and-shrink heuristic [" << timer << "]"
         << endl;
    cout << endl;
}

void MergeAndShrinkHeuristic::report_peak_memory_delta(bool final) const {
    if (final)
        cout << "Final";
    else
        cout << "Current";
    cout << " peak memory increase of merge-and-shrink computation: "
         << utils::get_peak_memory_in_kb() - starting_peak_memory << " KB"
         << endl;
}

void MergeAndShrinkHeuristic::dump_options() const {
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
}

void MergeAndShrinkHeuristic::warn_on_unusual_options() const {
    string dashes(79, '=');
    if (!label_reduction) {
        cerr << dashes << endl
             << "WARNING! You did not enable label reduction. This may "
            "drastically reduce the performance of merge-and-shrink!"
             << endl << dashes << endl;
    } else if (label_reduction->reduce_before_merging() && label_reduction->reduce_before_shrinking()) {
        cerr << dashes << endl
             << "WARNING! You set label reduction to be applied twice in "
            "each merge-and-shrink iteration, both before shrinking and\n"
            "merging. This double computation effort does not pay off "
            "for most configurations!"
             << endl << dashes << endl;
    } else {
        if (label_reduction->reduce_before_shrinking() &&
            (shrink_strategy->get_name() == "f-preserving"
             || shrink_strategy->get_name() == "random")) {
            cerr << dashes << endl
                 << "WARNING! Bucket-based shrink strategies such as "
                "f-preserving random perform best if used with label\n"
                "reduction before merging, not before shrinking!"
                 << endl << dashes << endl;
        }
        if (label_reduction->reduce_before_merging() &&
            shrink_strategy->get_name() == "bisimulation") {
            cerr << dashes << endl
                 << "WARNING! Shrinking based on bisimulation performs best "
                "if used with label reduction before shrinking, not\n"
                "before merging!"
                 << endl << dashes << endl;
        }
    }
}

pair<int, int> MergeAndShrinkHeuristic::compute_shrink_sizes(
    int size1, int size2) const {
    // Bound both sizes by max allowed size before merge.
    int new_size1 = min(size1, max_states_before_merge);
    int new_size2 = min(size2, max_states_before_merge);

    if (!utils::is_product_within_limit(new_size1, new_size2, max_states)) {
        int balanced_size = int(sqrt(max_states));

        if (new_size1 <= balanced_size) {
            // Size of the first transition system is small enough. Use whatever
            // is left for the second transition system.
            new_size2 = max_states / new_size1;
        } else if (new_size2 <= balanced_size) {
            // Inverted case as before.
            new_size1 = max_states / new_size2;
        } else {
            // Both transition systems are too big. We set both target sizes
            // to balanced_size. An alternative would be to set one to
            // N1 = balanced_size and the other to N2 = max_states /
            // balanced_size, to get closer to the allowed maximum.
            // However, this would make little difference (N2 would
            // always be N1, N1 + 1 or N1 + 2), and our solution has the
            // advantage of treating the transition systems symmetrically.
            new_size1 = balanced_size;
            new_size2 = balanced_size;
        }
    }
    assert(new_size1 <= size1 && new_size2 <= size2);
    assert(new_size1 <= max_states_before_merge);
    assert(new_size2 <= max_states_before_merge);
    assert(new_size1 * new_size2 <= max_states);
    return make_pair(new_size1, new_size2);
}

bool MergeAndShrinkHeuristic::shrink_transition_system(
    int index, int new_size) {
    assert(fts);
    const TransitionSystem &ts = fts->get_ts(index);
    assert(ts.is_solvable());
    int num_states = ts.get_size();
    if (num_states > min(new_size, shrink_threshold_before_merge)) {
        cout << ts.tag() << "current size: " << num_states;
        if (new_size < num_states)
            cout << " (new size limit: " << new_size;
        else
            cout << " (shrink threshold: " << shrink_threshold_before_merge;
        cout << ")" << endl;
        return shrink_strategy->shrink(*fts, index, new_size);
    }
    return false;
}

pair<bool, bool> MergeAndShrinkHeuristic::shrink_before_merge(
    int index1, int index2) {
    assert(fts);
    const TransitionSystem &ts1 = fts->get_ts(index1);
    const TransitionSystem &ts2 = fts->get_ts(index2);
    /*
      Compute the size limit for both transition systems as imposed by
      max_states and max_states_before_merge.
    */
    pair<int, int> new_sizes = compute_shrink_sizes(
        ts1.get_size(), ts2.get_size());

    /*
      For both transition systems, possibly compute and apply an
      abstraction.
      TODO: we could better use the given limit by increasing the size limit
      for the second shrinking if the first shrinking was larger than
      required.
    */
    bool shrunk1 = shrink_transition_system(index1, new_sizes.first);
    bool shrunk2 = shrink_transition_system(index2, new_sizes.second);
    return make_pair(shrunk1, shrunk2);
}

void MergeAndShrinkHeuristic::build_transition_system(const utils::Timer &timer) {
    // TODO: We're leaking memory here in various ways. Fix this.
    //       Don't forget that build_atomic_transition_systems also
    //       allocates memory.

    fts = utils::make_unique_ptr<FactoredTransitionSystem>(
        create_factored_transition_system(task_proxy, true, debug_transition_systems));
    print_time(timer, "after computation of atomic transition systems");
    cout << endl;
    int maximum_intermediate_size = 0;
    for (int i = 0; i < fts->get_size(); ++i) {
        int size = fts->get_ts(i).get_size();
        if (size > maximum_intermediate_size) {
            maximum_intermediate_size = size;
        }
    }

    unique_ptr<MergeStrategy> merge_strategy =
        merge_strategy_factory->compute_merge_strategy(task, *fts);

    vector<int> init_hvalue_increase;
    vector<int> remaining_labels;
    remaining_labels.push_back(fts->get_labels().compute_number_active_labels());
    int iteration_counter = 0;
    bool still_perfect = true;
    vector<pair<int, int>> merge_order;
    int final_index = -1; // TODO: get rid of this
    if (fts->is_solvable()) { // All atomic transition system are solvable.
        int number_of_merges = task_proxy.get_variables().size() - 1;
        for (int i = 0; i < number_of_merges; ++i) {
            // Choose next transition systems to merge
            pair<int, int> merge_indices = merge_strategy->get_next();
            int merge_index1 = merge_indices.first;
            int merge_index2 = merge_indices.second;
            cout << "Next pair of indices: (" << merge_index1 << ", " << merge_index2 << ")" << endl;
            assert(merge_index1 != merge_index2);
            merge_order.push_back(merge_indices);
            fts->statistics(merge_index1);
            fts->statistics(merge_index2);
            print_time(timer, "after computation of next merge");

            // Label reduction (before shrinking)
            if (label_reduction && label_reduction->reduce_before_shrinking()) {
                bool reduced = label_reduction->reduce(merge_indices, *fts);
                if (reduced) {
                    print_time(timer, "after label reduction");
                }
                remaining_labels.push_back(fts->get_labels().compute_number_active_labels());
            }

            // Shrinking
            pair<bool, bool> shrunk = shrink_before_merge(
                merge_index1, merge_index2);
            if (shrunk.first) {
                fts->statistics(merge_index1);
            }
            if (shrunk.second) {
                fts->statistics(merge_index2);
            }
            if (shrunk.first || shrunk.second) {
                print_time(timer, "after shrinking");
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

            // Label reduction (before merging)
            if (label_reduction && label_reduction->reduce_before_merging()) {
                bool reduced = label_reduction->reduce(merge_indices, *fts);
                if (reduced) {
                    print_time(timer, "after label reduction");
                }
                remaining_labels.push_back(fts->get_labels().compute_number_active_labels());
            }

            int init_dist1 = fts->get_init_state_goal_distance(merge_index1);
            int init_dist2 = fts->get_init_state_goal_distance(merge_index2);

            // Merging
            final_index = fts->merge(merge_index1, merge_index2);
            /*
              NOTE: both the shrinking strategy classes and the construction of
              the composite require input transition systems to be solvable.
            */
            if (!fts->is_solvable()) {
                break;
            }

            int abs_size = fts->get_ts(final_index).get_size();
            if (abs_size > maximum_intermediate_size) {
                maximum_intermediate_size = abs_size;
            }

            fts->statistics(final_index);
            print_time(timer, "after merging");

            int new_init_dist = fts->get_init_state_goal_distance(final_index);
            int difference = new_init_dist - max(init_dist1, init_dist2);
            cout << "Difference of init h values: " << difference << endl;
            init_hvalue_increase.push_back(difference);

            report_peak_memory_delta();
            cout << endl;
            ++iteration_counter;
        }
    }

    if (fts->is_solvable()) {
        cout << "Final transition system size: "
             << fts->get_ts(final_index).get_size() << endl;
        // need to finalize before calling "get_cost"
        fts->finalize();
        // TODO: after adopting the task interface everywhere, change this
        // back to compute_heuristic(task_proxy.get_initial_state())
        cout << "initial h value: "
             << fts->get_cost(task_proxy.get_initial_state())
             << endl;
    } else {
        cout << "Abstract problem is unsolvable!" << endl;
    }

    cout << "Maximum intermediate abstraction size: "
         << maximum_intermediate_size << endl;
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
    const vector<double> &pruning_statistics = fts->get_pruning_statistics();
    cout << "Relative pruning per iteration: " << pruning_statistics << endl;
    double summed_pruning = 0;
    for (double pruning : pruning_statistics) {
        summed_pruning += pruning;
    }
    // If pruning statistics are empty, then because the instance is unsolvable.
    // In this case, we return 0, which is the worst value possible for pruning.
    double average_pruning = 0;
    if (!pruning_statistics.empty()) {
        average_pruning =  summed_pruning / static_cast<double>(pruning_statistics.size());
    }
    cout << "Average relative pruning: " << average_pruning << endl;
    cout << "Iterations with merge tiebreaking: "
         << merge_strategy->get_iterations_with_tiebreaking() << endl;
    cout << "Total tiebreaking merge candidates: "
         << merge_strategy->get_total_tiebreaking_pair_count() << endl;

    merge_strategy = nullptr;
    shrink_strategy = nullptr;
    label_reduction = nullptr;
}

void MergeAndShrinkHeuristic::initialize() {
}

int MergeAndShrinkHeuristic::compute_heuristic(const GlobalState &global_state) {
    State state = convert_global_state(global_state);
    int cost = fts->get_cost(state);
    if (cost == -1)
        return DEAD_END;
    return cost;
}

void MergeAndShrinkHeuristic::add_shrink_limit_options_to_parser(OptionParser &parser) {
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

void MergeAndShrinkHeuristic::handle_shrink_limit_options_defaults(Options &opts) {
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
        cerr << "warning: max_states_before_merge exceeds max_states, "
             << "correcting." << endl;
        max_states_before_merge = max_states;
    }

    if (max_states < 1) {
        cerr << "error: transition system size must be at least 1" << endl;
        utils::exit_with(ExitCode::INPUT_ERROR);
    }

    if (max_states_before_merge < 1) {
        cerr << "error: transition system size before merge must be at least 1"
             << endl;
        utils::exit_with(ExitCode::INPUT_ERROR);
    }

    if (threshold == -1) {
        threshold = max_states;
    }
    if (threshold < 1) {
        cerr << "error: threshold must be at least 1" << endl;
        utils::exit_with(ExitCode::INPUT_ERROR);
    }
    if (threshold > max_states) {
        cerr << "warning: threshold exceeds max_states, correcting" << endl;
        threshold = max_states;
    }

    opts.set<int>("max_states", max_states);
    opts.set<int>("max_states_before_merge", max_states_before_merge);
    opts.set<int>("threshold_before_merge", threshold);
}

static Heuristic *_parse(OptionParser &parser) {
    parser.document_synopsis(
        "Merge-and-shrink heuristic",
        "This heuristic implements the algorithm described in the following "
        "paper:" + utils::format_paper_reference(
            {"Silvan Sievers", "Martin Wehrle", "Malte Helmert"},
            "Generalized Label Reduction for Merge-and-Shrink Heuristics",
            "http://ai.cs.unibas.ch/papers/sievers-et-al-aaai2014.pdf",
            "Proceedings of the 28th AAAI Conference on Artificial"
            " Intelligence (AAAI 2014)",
            "2358-2366",
            "AAAI Press 2014") + "\n" +
        "For a more exhaustive description of merge-and-shrink, see the journal "
        "paper" + utils::format_paper_reference(
            {"Malte Helmert", "Patrik Haslum", "Joerg Hoffmann", "Raz Nissim"},
            "Merge-and-Shrink Abstraction: A Method for Generating Lower Bounds"
            " in Factored State Spaces",
            "http://ai.cs.unibas.ch/papers/helmert-et-al-jacm2014.pdf",
            "Journal of the ACM 61 (3)",
            "16:1-63",
            "2014") + "\n" +
        "Please note that the journal paper describes the \"old\" theory of "
        "label reduction, which has been superseded by the above conference "
        "paper and is no longer implemented in Fast Downward.");
    parser.document_language_support("action costs", "supported");
    parser.document_language_support("conditional effects", "supported (but see note)");
    parser.document_language_support("axioms", "not supported");
    parser.document_property("admissible", "yes");
    parser.document_property("consistent", "yes");
    parser.document_property("safe", "yes");
    parser.document_property("preferred operators", "no");
    parser.document_note(
        "Note",
        "Conditional effects are supported directly. Note, however, that "
        "for tasks that are not factored (in the sense of the JACM 2014 "
        "merge-and-shrink paper), the atomic transition systems on which "
        "merge-and-shrink heuristics are based are nondeterministic, "
        "which can lead to poor heuristics even when only perfect shrinking "
        "is performed.");
    parser.document_note(
        "Note",
        "A currently recommended good configuration uses bisimulation "
        "based shrinking (selecting max states from 50000 to 200000 is "
        "reasonable), DFP merging, and the appropriate label "
        "reduction setting:\n"
        "merge_and_shrink(shrink_strategy=shrink_bisimulation(greedy=false),"
        "merge_strategy=merge_dfp(),label_reduction=label_reduction("
        "before_shrinking=true, before_merging=false),max_states=100000,"
        "threshold_before_merge=1)");

    // Merge strategy option.
    parser.add_option<shared_ptr<MergeStrategyFactory>>(
        "merge_strategy",
        "See detailed documentation for merge strategies. "
        "We currently recommend merge_dfp.");

    // Shrink strategy option.
    parser.add_option<shared_ptr<ShrinkStrategy>>(
        "shrink_strategy",
        "See detailed documentation for shrink strategies. "
        "We currently recommend shrink_bisimulation.");

    // Label reduction option.
    parser.add_option<shared_ptr<LabelReduction>>(
        "label_reduction",
        "See detailed documentation for labels. There is currently only "
        "one 'option' to use label_reduction. Also note the interaction "
        "with shrink strategies.",
        OptionParser::NONE);

    MergeAndShrinkHeuristic::add_shrink_limit_options_to_parser(parser);

    parser.add_option<bool>("debug_transition_systems", "store additional information "
                            "in transition systems for debug output.", "false");

    Heuristic::add_options_to_parser(parser);
    Options opts = parser.parse();
    if (parser.help_mode()) {
        return nullptr;
    }

    MergeAndShrinkHeuristic::handle_shrink_limit_options_defaults(opts);

    if (parser.dry_run()) {
        return nullptr;
    } else {
        return new MergeAndShrinkHeuristic(opts);
    }
}

static Plugin<Heuristic> _plugin("merge_and_shrink", _parse);
}
