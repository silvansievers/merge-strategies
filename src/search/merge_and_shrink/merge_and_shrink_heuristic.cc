#include "merge_and_shrink_heuristic.h"

#include "distances.h"
#include "factor_scoring_functions.h"
#include "factored_transition_system.h"
#include "merge_and_shrink_algorithm.h"
#include "merge_and_shrink_representation.h"
#include "transition_system.h"
#include "types.h"

#include "../option_parser.h"
#include "../plugin.h"

#include "../task_utils/task_properties.h"

#include "../utils/markup.h"
#include "../utils/system.h"

#include <cassert>
#include <iostream>
#include <utility>

using namespace std;
using utils::ExitCode;

namespace merge_and_shrink {
MergeAndShrinkHeuristic::MergeAndShrinkHeuristic(const options::Options &opts)
    : Heuristic(opts),
      verbosity(static_cast<Verbosity>(opts.get_enum("verbosity"))),
      partial_mas_method(static_cast<PartialMASMethod>(opts.get_enum("partial_mas_method"))) {
    if (opts.contains("factor_scoring_functions")) {
        factor_scoring_functions = opts.get_list<shared_ptr<FactorScoringFunction>>(
            "factor_scoring_functions");
    }

    cout << "Initializing merge-and-shrink heuristic..." << endl;
    MergeAndShrinkAlgorithm algorithm(opts);
    FactoredTransitionSystem fts = algorithm.build_factored_transition_system(task_proxy);
    finalize(fts);
    cout << "Done initializing merge-and-shrink heuristic." << endl << endl;
}

vector<int> get_remaining_candidates(
    const vector<int> &merge_candidates,
    const vector<int> &scores) {
    assert(merge_candidates.size() == scores.size());
    int best_score = -1;
    for (int score : scores) {
        if (score > best_score) {
            best_score = score;
        }
    }

    vector<int> result;
    for (size_t i = 0; i < scores.size(); ++i) {
        if (scores[i] == best_score) {
            result.push_back(merge_candidates[i]);
        }
    }
    return result;
}

int MergeAndShrinkHeuristic::find_best_factor(
    const FactoredTransitionSystem &fts) const {
    vector<int> current_indices;
    for (int index : fts) {
        current_indices.push_back(index);
    }

    for (const shared_ptr<FactorScoringFunction> &fsf : factor_scoring_functions) {
        vector<int> scores = fsf->compute_scores(fts, current_indices);
        current_indices = get_remaining_candidates(current_indices, scores);
        if (current_indices.size() == 1) {
            break;
        }
    }

    if (current_indices.size() != 1) {
        cerr << "More than one factor remained. Did you forget including "
                "a unique tie-breaking factor scoring function?" << endl;
        utils::exit_with(utils::ExitCode::SEARCH_CRITICAL_ERROR);
    }

    return current_indices.front();
}

void MergeAndShrinkHeuristic::finalize_factor(
    FactoredTransitionSystem &fts, int index) {
    auto final_entry = fts.extract_factor(index);
    unique_ptr<MergeAndShrinkRepresentation> mas_representation = move(final_entry.first);
    unique_ptr<Distances> distances = move(final_entry.second);
    if (!distances->are_goal_distances_computed()) {
        const bool compute_init = false;
        const bool compute_goal = true;
        distances->compute_distances(compute_init, compute_goal, verbosity);
    }
    assert(distances->are_goal_distances_computed());
    mas_representation->set_distances(*distances);
    mas_representations.push_back(move(mas_representation));
}

void MergeAndShrinkHeuristic::finalize(FactoredTransitionSystem &fts) {
    /*
      TODO: This method has quite a bit of fiddling with aspects of
      transition systems and the merge-and-shrink representation (checking
      whether distances have been computed; computing them) that we would
      like to have at a lower level. See also the TODO in
      factored_transition_system.h on improving the interface of that class
      (and also related classes like TransitionSystem etc).
    */

    int active_factors_count = fts.get_num_active_entries();
    if (verbosity >= Verbosity::NORMAL) {
        cout << "Number of remaining factors: " << active_factors_count << endl;
    }

    // Check if there is an unsolvable factor. If so, use it as heuristic.
    for (int index : fts) {
        if (!fts.is_factor_solvable(index)) {
            finalize_factor(fts, index);
            if (verbosity >= Verbosity::NORMAL) {
                cout << fts.get_transition_system(index).tag()
                     << "use this unsolvable factor as heuristic."
                     << endl;
            }
            return;
        }
    }

    // Iterate over all remaining factors and extract them.
    for (int index : fts) {
        finalize_factor(fts, index);
    }
    if (verbosity >= Verbosity::NORMAL) {
        cout << "Use all factors in a maximum heuristic." << endl;
    }

    assert(static_cast<int>(mas_representations.size()) == active_factors_count);
}

int MergeAndShrinkHeuristic::compute_heuristic(const GlobalState &global_state) {
    State state = convert_global_state(global_state);
    int heuristic = 0;
    for (const unique_ptr<MergeAndShrinkRepresentation> &mas_representation : mas_representations) {
        int cost = mas_representation->get_value(state);
        if (cost == PRUNED_STATE || cost == INF) {
            // If state is unreachable or irrelevant, we encountered a dead end.
            return DEAD_END;
        }
        heuristic = max(heuristic, cost);
    }
    return heuristic;
}

static shared_ptr<Heuristic> _parse(options::OptionParser &parser) {
    parser.document_synopsis(
        "Merge-and-shrink heuristic",
        "This heuristic implements the algorithm described in the following "
        "paper:" + utils::format_paper_reference(
            {"Silvan Sievers", "Martin Wehrle", "Malte Helmert"},
            "Generalized Label Reduction for Merge-and-Shrink Heuristics",
            "https://ai.dmi.unibas.ch/papers/sievers-et-al-aaai2014.pdf",
            "Proceedings of the 28th AAAI Conference on Artificial"
            " Intelligence (AAAI 2014)",
            "2358-2366",
            "AAAI Press 2014") + "\n" +
        "For a more exhaustive description of merge-and-shrink, see the journal "
        "paper" + utils::format_paper_reference(
            {"Malte Helmert", "Patrik Haslum", "Joerg Hoffmann", "Raz Nissim"},
            "Merge-and-Shrink Abstraction: A Method for Generating Lower Bounds"
            " in Factored State Spaces",
            "https://ai.dmi.unibas.ch/papers/helmert-et-al-jacm2014.pdf",
            "Journal of the ACM 61 (3)",
            "16:1-63",
            "2014") + "\n" +
        "Please note that the journal paper describes the \"old\" theory of "
        "label reduction, which has been superseded by the above conference "
        "paper and is no longer implemented in Fast Downward.\n\n"
        "The following paper describes how to improve the DFP merge strategy "
        "with tie-breaking, and presents two new merge strategies (dyn-MIASM "
        "and SCC-DFP):" + utils::format_paper_reference(
            {"Silvan Sievers", "Martin Wehrle", "Malte Helmert"},
            "An Analysis of Merge Strategies for Merge-and-Shrink Heuristics",
            "https://ai.dmi.unibas.ch/papers/sievers-et-al-icaps2016.pdf",
            "Proceedings of the 26th International Conference on Automated "
            "Planning and Scheduling (ICAPS 2016)",
            "294-298",
            "AAAI Press 2016") + "\n" +
        "Details of the algorithms and the implementation are described in the "
        "paper" + utils::format_paper_reference(
            {"Silvan Sievers"},
            "Merge-and-Shrink Heuristics for Classical Planning: Efficient "
            "Implementation and Partial Abstractions",
            "https://ai.dmi.unibas.ch/papers/sievers-socs2018.pdf",
            "Proceedings of the 11th Annual Symposium on Combinatorial Search "
            "(SoCS 2018)",
            "90-98",
            "AAAI Press 2018")
        );
    parser.document_language_support("action costs", "supported");
    parser.document_language_support("conditional effects", "supported (but see note)");
    parser.document_language_support("axioms", "not supported");
    parser.document_property("admissible", "yes (but see note)");
    parser.document_property("consistent", "yes (but see note)");
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
        "When pruning unreachable states, admissibility and consistency is "
        "only guaranteed for reachable states and transitions between "
        "reachable states. While this does not impact regular A* search which "
        "will never encounter any unreachable state, it impacts techniques "
        "like symmetry-based pruning: a reachable state which is mapped to an "
        "unreachable symmetric state (which hence is pruned) would falsely be "
        "considered a dead-end and also be pruned, thus violating optimality "
        "of the search.");
    parser.document_note(
        "Note",
        "When using a time limit on the main loop of the merge-and-shrink "
        "algorithm, the heuristic will compute the maximum over all heuristics "
        "induced by the remaining factors if terminating the merge-and-shrink "
        "algorithm early. Exception: if there is an unsolvable factor, it will "
        "be used as the exclusive heuristic since the problem is unsolvable.");
    parser.document_note(
        "Note",
        "A currently recommended good configuration uses bisimulation "
        "based shrinking, the merge strategy SCC-DFP, and the appropriate "
        "label reduction setting (max_states has been altered to be between "
        "10000 and 200000 in the literature):\n"
        "{{{\nmerge_and_shrink(shrink_strategy=shrink_bisimulation(greedy=false),"
        "merge_strategy=merge_sccs(order_of_sccs=topological,merge_selector="
        "score_based_filtering(scoring_functions=[goal_relevance,dfp,"
        "total_order])),label_reduction=exact(before_shrinking=true,"
        "before_merging=false),max_states=50000,threshold_before_merge=1)\n}}}\n"
        "Note that for versions of Fast Downward prior to 2016-08-19, the "
        "syntax differs. See the recommendation in the file "
        "merge_and_shrink_heuristic.cc for an example configuration.");

    Heuristic::add_options_to_parser(parser);
    add_merge_and_shrink_algorithm_options_to_parser(parser);

    vector<string> partial_mas_method;
    vector<string> partial_mas_method_docs;
    partial_mas_method.push_back("none");
    partial_mas_method_docs.push_back(
        "none: attempt to compute a merge-and-shrink abstraction over all "
        "variables of the planning task. Do not set a finite value for any of"
        "the options main_loop_max_time or num_transitions_to_abort");
    partial_mas_method.push_back("single");
    partial_mas_method_docs.push_back(
        "single: choose a single factor of the remaining factors to serve as"
        "the abstraction for the heuristic. The factor is chosen according to "
        "the factor scoring functions provided via factor_scoring_functions.");
    partial_mas_method.push_back("maximum");
    partial_mas_method_docs.push_back(
        "maximum: retain all remaining factors and compute the maximum "
        "heuristic over all these abstractions.");
    parser.add_enum_option(
        "partial_mas_method",
        partial_mas_method,
        "Method to determine the final heuristic given an early abortion, "
        "such as due to reaching the time or transitions limit.",
        "none",
        partial_mas_method_docs);
    parser.add_list_option<shared_ptr<FactorScoringFunction>>(
        "factor_scoring_functions",
        "The list of factor scoring functions used to compute scores for "
        "remaining factors if computing partial merge-and-shrink abstractions, "
        "i.e., if partial_mas_method != none.",
        options::OptionParser::NONE);

    options::Options opts = parser.parse();
    if (parser.help_mode()) {
        return nullptr;
    }

    handle_shrink_limit_options_defaults(opts);

    if (parser.dry_run()) {
        double main_loop_max_time = opts.get<double>("main_loop_max_time");
        int num_transitions_to_abort = opts.get<int>("num_transitions_to_abort");
        PartialMASMethod partial_mas_method = static_cast<PartialMASMethod>(opts.get_enum("partial_mas_method"));
        if (partial_mas_method != PartialMASMethod::None
            && (main_loop_max_time == numeric_limits<double>::infinity()
                && num_transitions_to_abort == INF)) {
            cerr << "If using a partial merge-and-shrink method, you must "
                "use a finite value for at least one of main_loop_max_time and "
                "num_transitions_to_abort. "
                 << endl;
            utils::exit_with(utils::ExitCode::SEARCH_INPUT_ERROR);
        }
        if (partial_mas_method == PartialMASMethod::None
            && (main_loop_max_time < INF || num_transitions_to_abort < INF)) {
            cerr << "If using a finite value to any of main_loop_max_time and "
                "num_transitions_to_abort, you also must use a partial "
                "merge-and-shrink method."
                 << endl;
            utils::exit_with(utils::ExitCode::SEARCH_INPUT_ERROR);
        }
        if (partial_mas_method == PartialMASMethod::Single
            && !opts.contains("factor_scoring_functions")) {
            cerr << "If using the partial merge-and-shrink method single, "
                "you must specify a least one factor scoring function!"
                 << endl;
            utils::exit_with(utils::ExitCode::SEARCH_INPUT_ERROR);
        }
        return nullptr;
    } else {
        return make_shared<MergeAndShrinkHeuristic>(opts);
    }
}

static options::Plugin<Evaluator> _plugin("merge_and_shrink", _parse);
}
