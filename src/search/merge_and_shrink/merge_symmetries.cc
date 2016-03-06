#include "merge_symmetries.h"

#include "factored_transition_system.h"
#include "merge_linear.h"
#include "merge_dfp.h"
#include "merge_miasm.h"
#include "transition_system.h"

#include "miasm/merge_tree.h"

#include "symmetries/symmetry_group.h"

#include "../globals.h"
#include "../option_parser.h"
#include "../plugin.h"
#include "../variable_order_finder.h"

#include "../utils/logging.h"
#include "../utils/rng.h"

#include <algorithm>
#include <iomanip>

using namespace std;

namespace merge_and_shrink {
MergeSymmetries::MergeSymmetries(const Options &options_)
    : MergeStrategy(),
      max_bliss_iterations(options_.get<int>("max_bliss_iterations")),
      bliss_call_time_limit(options_.get<int>("bliss_call_time_limit")),
      bliss_remaining_time_budget(options_.get<int>("bliss_total_time_budget")),
      fallback_strategy(FallbackStrategy(options_.get_enum("fallback_strategy"))),
      iteration_counter(0),
      number_of_applied_symmetries(0),
      bliss_limit_reached(false),
      pure_fallback_strategy(true),
      merge_dfp(nullptr),
      miasm_merge_tree(nullptr) {
    options = new Options(options_);
}

MergeSymmetries::~MergeSymmetries() {
    delete merge_dfp;
    delete options;
    delete miasm_merge_tree;
}

void MergeSymmetries::dump_statistics() {
    if (!bliss_times.empty()) {
        sort(bliss_times.begin(), bliss_times.end());
        double summed_up_bliss_times = 0;
        size_t total_bliss_calls = bliss_times.size();
        cout << "Total bliss calls: " << total_bliss_calls << endl;
        for (size_t i = 0; i < total_bliss_calls; ++i) {
            summed_up_bliss_times += bliss_times[i];
        }
        double average_bliss_time = summed_up_bliss_times
                / (double) total_bliss_calls;
        cout << setprecision(5) << fixed << "Average bliss time: " << average_bliss_time << endl;
        double median_bliss_time;
        if (total_bliss_calls % 2 == 0) {
            size_t lower_median_index = (total_bliss_calls - 1) / 2;
            size_t higher_median_index = lower_median_index + 1;
            median_bliss_time = (bliss_times[lower_median_index]
                                 + bliss_times[higher_median_index])
                    / 2.0;
        } else {
            size_t median_index = total_bliss_calls / 2;
            median_bliss_time = bliss_times[median_index];
        }
        cout << setprecision(5) << fixed << "Median bliss time: " << median_bliss_time << endl;
    }
}

void MergeSymmetries::dump_strategy_specific_options() const {
    cout << "Options for merge symmetries:" << endl;
    cout << "    symmetries for shrinking: ";
    int symmetries_for_shrinking = options->get_enum("symmetries_for_shrinking");
    switch (symmetries_for_shrinking) {
        case 0: {
            cout << "none";
            break;
        }
        case 1: {
            cout << "atomic";
            break;
        }
        case 2: {
            cout << "local";
            break;
        }
    }
    cout << endl;
    cout << "    symmetries for merging: ";
    int symmetries_for_merging = options->get_enum("symmetries_for_merging");
    switch (symmetries_for_merging) {
        case 0: {
            cout << "none";
            break;
        }
        case 1: {
            cout << "smallest";
            break;
        }
        case 2: {
            cout << "largest";
            break;
        }
    }
    cout << endl;
    if (symmetries_for_merging) {
        cout << "    external merging: ";
        switch (options->get_enum("external_merging")) {
            case 0: {
                cout << "merge for atomic symmetry";
                break;
            }
            case 1: {
                cout << "merge for local symmetry";
                break;
            }
        }
        cout << endl;
        cout << "    internal merging: ";
        switch (options->get_enum("internal_merging")) {
            case 0: {
                cout << "linear";
                break;
            }
            case 1: {
                cout << "non linear";
                break;
            }
        }
        cout << endl;
    }
    cout << "    maxium number of m&s iterations with bliss: "
         << max_bliss_iterations << endl;
    cout << "    time limit for single bliss calls (0 means unlimited): "
         << bliss_call_time_limit << endl;
    cout << "    total time budget for bliss (0 means unlimited): "
         << options->get<int>("bliss_total_time_budget") << endl;
    cout << "    stop searching for symmetries once no symmetry was found: "
         << (options->get<bool>("stop_after_no_symmetries") ? "yes" : "no") << endl;
    cout << "    stabilize transition systems: "
         << (options->get<bool>("stabilize_transition_systems") ? "yes" : "no") << endl;
    cout << "    fallback merge strategy: ";
    switch (fallback_strategy) {
    case LINEAR:
        cout << "linear";
        break;
    case DFP:
        cout << "dfp";
        break;
    case MIASM:
        cout << "miasm";
        break;
    }
    cout << endl;
}

void MergeSymmetries::initialize(const shared_ptr<AbstractTask> task) {
    MergeStrategy::initialize(task);
    if (fallback_strategy == LINEAR) {
        TaskProxy task_proxy(*task);
        VariableOrderFinder order(task,
                                  VariableOrderType(options->get_enum("variable_order")));
        linear_merge_order.reserve(task_proxy.get_variables().size());
        while (!order.done()) {
            linear_merge_order.push_back(order.next());
        }
    } else if (fallback_strategy == DFP) {
        merge_dfp = new MergeDFP(*options);
        merge_dfp->initialize(task);
    } else if (fallback_strategy == MIASM) {
        MergeMiasm merge_miasm(*options);
        miasm_merge_tree = merge_miasm.compute_merge_tree(task);
    } else {
        ABORT("unknown fallback merge strategy");
    }
}

pair<int, int> MergeSymmetries::get_next_miasm(FactoredTransitionSystem &fts) {
    tree<set<int>> &merge_tree = miasm_merge_tree->get_tree();
    pair<int, int> next_indices;

    // TODO: original miasm uses sibling iterators. is this the same tree
    // traversal order in as using leaf iterators in our case?
    for (tree<set<int>>::leaf_iterator i_node = merge_tree.begin_leaf();
            i_node != merge_tree.end_leaf(); i_node++) {
        for (tree<set<int>>::leaf_iterator j_node = i_node;
                j_node != merge_tree.end_leaf(); j_node++) {
            if (merge_tree.parent(i_node) == merge_tree.parent(j_node)
                    && i_node != j_node) {
                set<int> &var_set1 = *i_node;
                set<int> &var_set2 = *j_node;
                assert(var_set1.size() == 1);
                assert(var_set2.size() == 1);
                int first_index = *var_set1.begin();
                int second_index = *var_set2.begin();
                next_indices.first = first_index;
                next_indices.second = second_index;

                int merged_index = fts.get_size();
                set<int> &parent_var_set = *merge_tree.parent(i_node);
                parent_var_set.clear();
                parent_var_set.insert(merged_index);
                merge_tree.erase(i_node);
                merge_tree.erase(j_node);
                return next_indices;
            }
        }
    }
    next_indices.first = -1;
    next_indices.second = -1;
    ABORT("Did not find MIASM next pair!");
    return next_indices;
}

void MergeSymmetries::update_miasm_merge_tree(
    FactoredTransitionSystem &fts,
    const pair<int, int> &next_merge) {
    tree<set<int>> &merge_tree = miasm_merge_tree->get_tree();
    // Upate merge tree:
//    cout << "tree before update" << endl;
//    cout << miasm_merge_tree << endl;

    // 1) Find the two indices in the tree
    int first_index = next_merge.first;
    int second_index = next_merge.second;
    tree<set<int>>::pre_order_iterator *first_index_it = 0;
    tree<set<int>>::pre_order_iterator *second_index_it = 0;
    for (tree<set<int>>::pre_order_iterator node_it = merge_tree.begin();
            node_it != merge_tree.end(); node_it++) {
        const set<int> &node = *node_it;
        if (node.count(first_index)) {
            first_index_it = new tree<set<int>>::pre_order_iterator(node_it);
        }
        if (node.count(second_index)) {
            second_index_it = new tree<set<int>>::pre_order_iterator(node_it);
        }
    }
    assert(first_index_it);
    assert(second_index_it);
    assert(first_index_it != second_index_it);

    int merged_index = fts.get_size();
    if (merge_tree.parent(*first_index_it) == merge_tree.parent(*second_index_it)) {
        // MIASM would merge the two indices anyway
        set<int> &parent_var_set = *merge_tree.parent(*first_index_it);
        parent_var_set.clear();
        parent_var_set.insert(merged_index);
        merge_tree.erase(*first_index_it);
        merge_tree.erase(*second_index_it);
    } else {
        // Merge the two tree nodes and arbitrarily assign the merged node to
        // one of the nodes' parent, and remove the other node's parent.
        tree<set<int>>::pre_order_iterator root = merge_tree.begin();
        assert(merge_tree.depth(root) == 0);
        assert(*first_index_it != root);
        assert(*second_index_it != root);

        // TODO: probably need to find existing node fts.get_size(); or is it in parent?

        tree<set<int>>::pre_order_iterator *chosen_survivor_node;
        if (*first_index_it == root) {
            chosen_survivor_node = first_index_it;
        } else if (*second_index_it == root) {
            chosen_survivor_node = second_index_it;
        } else if (merge_tree.depth(*first_index_it) < merge_tree.depth(*second_index_it)){
            chosen_survivor_node = first_index_it;
        } else if (merge_tree.depth(*first_index_it) > merge_tree.depth(*second_index_it)) {
            chosen_survivor_node = second_index_it;
        } else {
            int random = g_rng(2);
            cout << "random: " << random << endl;
            chosen_survivor_node = (random ? second_index_it : first_index_it);
        }
        tree<set<int>>::pre_order_iterator *to_be_erased_node =
                (chosen_survivor_node == first_index_it ? second_index_it : first_index_it);
//        cout << "chosen survivor: " << *((*chosen_survivor_node)->begin())
//             << " with depth " << tree.depth(*chosen_survivor_node) << endl;
//        cout << "chosen to be deleted: " << *((*to_be_erased_node)->begin())
//             << " with depth " << tree.depth(*to_be_erased_node) << endl;

        assert(merge_tree.depth(*to_be_erased_node) >= 2);

        // update survivor node
        (*chosen_survivor_node)->clear();
        (*chosen_survivor_node)->insert(fts.get_size());

        tree<set<int>>::pre_order_iterator to_be_erased_nodes_parent = merge_tree.parent(*to_be_erased_node);
        tree<set<int>>::pre_order_iterator to_be_erased_nodes_parents_parent = merge_tree.parent(to_be_erased_nodes_parent);
        merge_tree.erase((*to_be_erased_node));
        merge_tree.reparent(to_be_erased_nodes_parents_parent, to_be_erased_nodes_parent);
        merge_tree.erase(to_be_erased_nodes_parent);
    }

    delete first_index_it;
    delete second_index_it;

//    cout << "tree after update" << endl;
//    cout << miasm_merge_tree << endl;
}

pair<int, int> MergeSymmetries::get_next(FactoredTransitionSystem &fts) {
    assert(!done());
    ++iteration_counter;

    if (!bliss_limit_reached && iteration_counter <= max_bliss_iterations && merge_order.empty()) {
        double time_limit = 0.0;
        if (bliss_remaining_time_budget) {
            // we round everything down to the full second
            time_limit = bliss_remaining_time_budget;
        } else if (bliss_call_time_limit) {
            time_limit = bliss_call_time_limit;
        }
        cout << "Setting bliss time limit to " << time_limit << endl;
        options->set<double>("bliss_time_limit", time_limit);
        SymmetryGroup symmetry_group(*options);
        bool applied_symmetries =
            symmetry_group.find_and_apply_symmetries(fts, merge_order);
        if (applied_symmetries) {
            ++number_of_applied_symmetries;
        }
        if (symmetry_group.is_bliss_limit_reached()) {
            bliss_limit_reached = true;
        }
        if (pure_fallback_strategy && (!merge_order.empty())) {
            pure_fallback_strategy = false;
            cout << "not pure fallback strategy anymore" << endl;
        }
        double bliss_time = symmetry_group.get_bliss_time();
        bliss_times.push_back(bliss_time);
        if (bliss_remaining_time_budget) {
            bliss_remaining_time_budget -= bliss_time;
            if (bliss_remaining_time_budget <= 0) {
                assert(bliss_limit_reached);
                // in case of measurement inaccuracies, set bliss limit
                // reached to true in any case.
                bliss_limit_reached = true;
            }
            cout << "Remaining bliss time budget " << bliss_remaining_time_budget << endl;
        }
        cout << "Number of applied symmetries: " << number_of_applied_symmetries << endl;
    }

    if (remaining_merges == 1) {
        dump_statistics();
    }

    if (merge_order.empty()) {
        if (fallback_strategy == LINEAR) {
            // Linear Strategy
            int first = linear_merge_order[0];
            linear_merge_order.erase(linear_merge_order.begin());
            int second = linear_merge_order[0];
            linear_merge_order[0] = fts.get_size();
            cout << "Next pair (linear strategy): " << first << ", " << second << endl;
            --remaining_merges;
            return make_pair(first, second);
        } else if (fallback_strategy == DFP) {
            --remaining_merges;
            return merge_dfp->get_next(fts);
        } else if (fallback_strategy == MIASM) {
            --remaining_merges;
            return get_next_miasm(fts);
        } else {
            ABORT("unknown fallback merge strategy");
        }
    }

    assert(false);

    pair<int, int> next_merge = merge_order.front();
    if (!fts.is_active(next_merge.first) || !fts.is_active(next_merge.second)) {
        cerr << "Problem with the merge strategy: invalid indices" << endl;
        utils::exit_with(utils::ExitCode::CRITICAL_ERROR);
    }
    merge_order.erase(merge_order.begin());

    if (fallback_strategy == LINEAR) {
        // Update linear merge order in case we need to fallback on it
        bool found_entry = false;
        size_t counter = 0;
        for (vector<int>::iterator it = linear_merge_order.begin();
             it != linear_merge_order.end(); ++it) {
            if (!found_entry && (*it == next_merge.first ||
                                 *it == next_merge.second)) {
                //cout << "updating entry " << counter << " (" << *it << ")to " << fts.get_size(); << endl;
                found_entry = true;
                linear_merge_order[counter] = fts.get_size();
            } else if (found_entry && (*it == next_merge.first ||
                                       *it == next_merge.second)) {
                //cout << "erasing entry " << counter << " (" << *it << ")" << endl;
                linear_merge_order.erase(it);
                break;
            }
            ++counter;
        }
        assert(found_entry);
    } else if (fallback_strategy == MIASM) {
        // Update miasm merge order for future fallback cases
        update_miasm_merge_tree(fts, next_merge);
    }

    --remaining_merges;
    return next_merge;
}

string MergeSymmetries::name() const {
    return "symmetries";
}

static shared_ptr<MergeStrategy> _parse(OptionParser &parser) {
    // Options for symmetries computation
    parser.add_option<int>("max_bliss_iterations", "maximum ms iteration until "
                           "which bliss is allowed to run.",
                           "infinity");
    parser.add_option<int>("bliss_call_time_limit", "time in seconds one bliss "
                           "run is allowed to last at most (0 means no limit)",
                           "0");
    parser.add_option<int>("bliss_total_time_budget", "time in seconds bliss is "
                           "allowed to run overall (0 means no limit)",
                           "0");
    parser.add_option<bool>("stop_after_no_symmetries", "stop calling bliss "
                            "after unsuccesfull previous bliss call.",
                           "False");
    vector<string> symmetries_for_shrinking;
    symmetries_for_shrinking.push_back("NO_SHRINKING");
    symmetries_for_shrinking.push_back("ATOMIC");
    symmetries_for_shrinking.push_back("LOCAL");
    parser.add_enum_option("symmetries_for_shrinking",
                           symmetries_for_shrinking,
                           "choose the type of symmetries used for shrinking: "
                           "no shrinking, "
                           "only atomic symmetries, "
                           "local symmetries.",
                           "NO_SHRINKING");
    vector<string> symmetries_for_merging;
    symmetries_for_merging.push_back("NO_MERGING");
    symmetries_for_merging.push_back("SMALLEST");
    symmetries_for_merging.push_back("LARGEST");
    parser.add_enum_option("symmetries_for_merging",
                           symmetries_for_merging,
                           "choose the type of symmetries that should determine "
                           "the set of transition systems to be merged: "
                           "the smallest or the largest",
                           "SMALLEST");
    vector<string> external_merging;
    external_merging.push_back("MERGE_FOR_ATOMIC");
    external_merging.push_back("MERGE_FOR_LOCAL");
    parser.add_enum_option("external_merging",
                           external_merging,
                           "choose the set of transition systems to be merged: "
                           "merge for atomic: merge all transition systems affected "
                           "by the chosen symmetry, or "
                           "merge for local: merge only the transition systems "
                           "mapped (in cycles) to others. only merge every "
                           "cycle separately.",
                           "MERGE_FOR_ATOMIC");
    vector<string> internal_merging;
    internal_merging.push_back("LINEAR");
    internal_merging.push_back("NON_LINEAR");
    parser.add_enum_option("internal_merging",
                           internal_merging,
                           "choose the order in which to merge the set of "
                           "transition systems to be merged (only useful with "
                           "MERGE_FOR_ATOMIC): "
                           "linear (obvious), "
                           "non linear, which means to first merge every cycle, "
                           "and then the resulting intermediate transition systems.",
                           "LINEAR");

    // Options for GraphCreator
    parser.add_option<bool>("stabilize_transition_systems", "compute symmetries that "
                            "stabilize transition systems, i.e. that are local.", "false");
    parser.add_option<bool>("debug_graph_creator", "produce dot readable output "
                            "from the graph generating methods", "false");

    // Options for fallback merge strategy
    vector<string> fallback_strategy;
    fallback_strategy.push_back("linear");
    fallback_strategy.push_back("dfp");
    fallback_strategy.push_back("miasm");
    parser.add_enum_option("fallback_strategy",
                           fallback_strategy,
                           "choose a merge strategy: linear (specify "
                           "variable_order), dfp, or miasm.",
                           "dfp");

    MergeLinear::add_options_to_parser(parser);
    MergeDFP::add_options_to_parser(parser);
    MergeMiasm::add_options_to_parser(parser);

    Options options = parser.parse();
    if (options.get<int>("bliss_call_time_limit")
            && options.get<int>("bliss_total_time_budget")) {
        cerr << "Please only specify bliss_call_time_limit or "
                "bliss_total_time_budget but not both" << endl;
        utils::exit_with(utils::ExitCode::CRITICAL_ERROR);
    }
    if (options.get_enum("symmetries_for_shrinking") == 0
            && options.get_enum("symmetries_for_merging") == 0) {
        cerr << "Please use symmetries at least for shrinking or merging." << endl;
        utils::exit_with(utils::ExitCode::CRITICAL_ERROR);
    }
    if (parser.dry_run())
        return nullptr;
    else
        return make_shared<MergeSymmetries>(options);
}

static PluginShared<MergeStrategy> _plugin("merge_symmetries", _parse);
}
