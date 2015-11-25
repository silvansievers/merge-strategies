#include "merge_sccs.h"

#include "factored_transition_system.h"
#include "merge_dfp.h"
#include "transition_system.h"

#include "../causal_graph.h"
#include "../option_parser.h"
#include "../plugin.h"
#include "../scc.h"
#include "../variable_order_finder.h"

#include <algorithm>
#include <cassert>
#include <iostream>

using namespace std;

bool compare_sccs_increasing(const vector<int> &lhs, const vector<int> &rhs) {
    return lhs.size() < rhs.size();
}

bool compare_sccs_decreasing(const vector<int> &lhs, const vector<int> &rhs) {
    return lhs.size() > rhs.size();
}

MergeSCCs::MergeSCCs(const Options &options_)
    : MergeStrategy(),
      order_of_sccs(OrderOfSCCs(options_.get_enum("order_of_sccs"))),
      internal_merge_order(InternalMergeOrder(options_.get_enum("internal_merge_order"))),
      merged_sccs_merge_order(MergedSCCsMergeOrder(options_.get_enum("merged_sccs_merge_order"))),
      merge_dfp(nullptr),
      number_of_merges_for_scc(0),
      merged_all_sccs(false),
      start_merging_sccs(true) {
    options = new Options(options_);
}

MergeSCCs::~MergeSCCs() {
    delete merge_dfp;
}

void MergeSCCs::initialize(const std::shared_ptr<AbstractTask> task) {
    MergeStrategy::initialize(task);

    TaskProxy task_proxy(*task);
    VariablesProxy vars = task_proxy.get_variables();
    int num_vars = vars.size();

    if (internal_merge_order == DFP1 || merged_sccs_merge_order == DFP2) {
        merge_dfp = new MergeDFP(*options);
        merge_dfp->initialize(task);
    }
    if (internal_merge_order == LINEAR1 || merged_sccs_merge_order == LINEAR2) {
        VariableOrderFinder vof(task, VariableOrderType(options->get_enum("variable_order")));
        linear_variable_order.reserve(num_vars);
        while (!vof.done()) {
            linear_variable_order.push_back(vof.next());
        }
        cout << "linear variable order: " << linear_variable_order << endl;
    }
    delete options;
    options = nullptr;

    // Compute SCCs of the causal graph
    vector<vector<int>> cg;
    cg.reserve(num_vars);
    if (internal_merge_order == LINEAR_MANUAL) {
        is_causally_linked.resize(num_vars, false);
        is_goal_variable.resize(num_vars, false);
        for (FactProxy goal : task_proxy.get_goals())
            is_goal_variable[goal.get_variable().get_id()] = true;
        cg_successors.reserve(num_vars);
        cg_predecessors.reserve(num_vars);
    }
    for (VariableProxy var : vars) {
        const std::vector<int> &successors =
            task_proxy.get_causal_graph().get_successors(var.get_id());
        cg.push_back(successors);
        if (internal_merge_order == LINEAR_MANUAL) {
            cg_successors.push_back(successors);
            cg_predecessors.push_back(task_proxy.get_causal_graph().get_predecessors(var.get_id()));
        }
    }
//    cout << "CG:" << endl;
//    for (size_t var = 0; var < cg.size(); ++var) {
//        cout << var << ": " << cg[var] << endl;
//    }
    SCC scc(cg);
    vector<vector<int>> sccs(scc.get_result());

    // Put the SCCs in the desired order
    switch (order_of_sccs) {
    case TOPOLOGICAL:
        // SCCs are computed in topological order
        break;
    case REVERSE_TOPOLOGICAL:
        // SCCs are computed in topological order
        reverse(sccs.begin(), sccs.end());
        break;
    case DECREASING:
        sort(sccs.begin(), sccs.end(), compare_sccs_decreasing);
        break;
    case INCREASING:
        sort(sccs.begin(), sccs.end(), compare_sccs_increasing);
        break;
    }

    /*
      Compute the indices at which the merged SCCs can be found when all
      SCCs have been merged.
    */
    int index = num_vars - 1;
    cout << "found cg sccs:" << endl;
    indices_of_merged_sccs.reserve(sccs.size());
    int largest_scc_size = 0;
    for (const vector<int> &scc : sccs) {
        cout << scc << endl;
        int scc_size = scc.size();
        if (scc_size > largest_scc_size) {
            largest_scc_size = scc_size;
        }
        if (scc_size == 1) {
            indices_of_merged_sccs.push_back(scc.front());
        } else {
            index += scc_size - 1;
            indices_of_merged_sccs.push_back(index);
            // only store non-singleton sccs for internal merging
            cg_sccs.push_back(unordered_set<int>(scc.begin(), scc.end()));
        }
    }
    if (sccs.size() == 1) {
        cout << "Only one single SCC" << endl;
    }
    if (static_cast<int>(sccs.size()) == num_vars) {
        cout << "Only singleton SCCs" << endl;
        assert(cg_sccs.empty());
        merged_all_sccs = true;
    }
    cout << "indices of merged sccs: " << indices_of_merged_sccs << endl;
    current_scc_ts_indices.reserve(largest_scc_size);
}

pair<int, int> MergeSCCs::get_next_linear(
    const shared_ptr<FactoredTransitionSystem> fts,
    const vector<int> available_indices,
    int most_recent_index,
    bool two_indices) const {
    int next_index1 = -1;
    if (!two_indices) {
        next_index1 = most_recent_index;
    }
    int next_index2 = -1;
    for (int var : linear_variable_order) {
        for (int index : available_indices) {
            if (index != next_index1) {
                const vector<int> &incorporated_variables =
                    fts->get_ts(index).get_incorporated_variables();
                vector<int>::const_iterator it = find(incorporated_variables.begin(),
                                                incorporated_variables.end(),
                                                var);
                if (it != incorporated_variables.end()) { // ts contains var
                    if (next_index1 == -1) {
                        next_index1 = index;
                        break;
                    } else {
                        assert(next_index2 == -1);
                        next_index2 = index;
                        break;
                    }
                }
            }
        }
        if (next_index1 != -1 && next_index2 != -1) {
            break;
        }
    }
    return make_pair(next_index1, next_index2);
}

pair<int, int> MergeSCCs::get_next(
    shared_ptr<FactoredTransitionSystem> fts) {
    assert(!done());

    pair<int, int > next_pair = make_pair(-1, -1);
    int most_recent_index = fts->get_size() - 1;
    if (!merged_all_sccs) {

        bool first_merge_of_scc = false; // needed for linear merging
        if (!number_of_merges_for_scc) {
            // We did not start merging a specific SCC yet
            first_merge_of_scc = true;
            assert(!cg_sccs.empty());
            const unordered_set<int> &current_scc = cg_sccs.front();
            assert(current_scc.size() > 1);
            assert(current_scc_ts_indices.empty());
            number_of_merges_for_scc = current_scc.size() - 1;
            // Initialize current transition systems with all those contained in the scc
            current_scc_ts_indices.insert(current_scc_ts_indices.end(),
                                          current_scc.begin(), current_scc.end());
        } else {
            // Add the newest transition system to the set of current ones of the scc
            current_scc_ts_indices.push_back(most_recent_index);
        }

        if (number_of_merges_for_scc > 1) {
            if (internal_merge_order == LINEAR1) {
                next_pair = get_next_linear(fts,
                                            current_scc_ts_indices,
                                            most_recent_index,
                                            first_merge_of_scc);
            } else if (internal_merge_order == DFP1) {
                next_pair = merge_dfp->get_next(fts, current_scc_ts_indices);
            } else if (internal_merge_order == LINEAR_MANUAL) {
                int next_index1 = -1;
                int next_index2 = -1;
                if (first_merge_of_scc) {
                    int max_num_vars = is_causally_linked.size();
                    is_causally_linked.assign(max_num_vars, false);
                    for (int var : current_scc_ts_indices) {
                        if (is_goal_variable[var]) {
                            next_index1 = var;
                            break;
                        }
                    }
                    if (next_index1 == -1) {
                        next_index1 = current_scc_ts_indices.front();
                    }
                    const vector<int> &successors = cg_successors[next_index1];
                    for (int succ : successors) {
                        is_causally_linked[succ] = true;
                    }
                    const vector<int> &predecessors = cg_predecessors[next_index1];
                    for (int pred : predecessors) {
                        is_causally_linked[pred] = true;
                    }
                } else {
                    next_index1 = most_recent_index;
                }
                assert(next_index1 != -1);
                for (int var : current_scc_ts_indices) {
                    if (var != next_index1 && is_causally_linked[var]) {
                        next_index2 = var;
                        break;
                    }
                }
                assert(next_index2 != -1);
                const vector<int> &successors = cg_successors[next_index2];
                for (int succ : successors) {
                    is_causally_linked[succ] = true;
                }
                const vector<int> &predecessors = cg_predecessors[next_index2];
                for (int pred : predecessors) {
                    is_causally_linked[pred] = true;
                }
                next_pair = make_pair(next_index1, next_index2);
            }

            // Remove the two merged indices from the current set of indices
            for (vector<int>::iterator it = current_scc_ts_indices.begin();
                it != current_scc_ts_indices.end(); ) {
                if (*it == next_pair.first || *it == next_pair.second) {
                    it = current_scc_ts_indices.erase(it);
                } else {
                    ++it;
                }
            }
        } else {
            assert(number_of_merges_for_scc == 1);
            assert(current_scc_ts_indices.size() == 2);
            next_pair = make_pair(current_scc_ts_indices[0],
                current_scc_ts_indices[1]);

            current_scc_ts_indices.clear();
            cg_sccs.erase(cg_sccs.begin());
            if (cg_sccs.empty()) {
                merged_all_sccs = true;
            }
        }

        --number_of_merges_for_scc;
    } else {
        // TODO: now this is very similar to the phase of internally merging sccs.
        // need to reduce code duplication

        bool first_merge_of_merged_sccs_merging = false; // needed for linear merging
        if (start_merging_sccs) {
            // If we end up here the first time, we have at least 2 SCCs
            assert(indices_of_merged_sccs.size() > 1);
            start_merging_sccs = false;
            first_merge_of_merged_sccs_merging = true;
        }  else {
            // Add the newest transition system to the set of current merged
            // sccs indices
            indices_of_merged_sccs.push_back(most_recent_index);
        }

        if (indices_of_merged_sccs.size() == 2) {
            next_pair = make_pair(indices_of_merged_sccs[0],
                indices_of_merged_sccs[1]);
            indices_of_merged_sccs.clear();
        } else {
            if (merged_sccs_merge_order == LINEAR2) {
                next_pair = get_next_linear(fts,
                                            indices_of_merged_sccs,
                                            most_recent_index,
                                            first_merge_of_merged_sccs_merging);
            } else if (merged_sccs_merge_order == DFP2) {
                next_pair = merge_dfp->get_next(fts, indices_of_merged_sccs);
            }

            // Remove the two indices from indices_of_merged
            for (vector<int>::iterator it = indices_of_merged_sccs.begin();
                it != indices_of_merged_sccs.end(); ) {
                if (*it == next_pair.first || *it == next_pair.second) {
                    it = indices_of_merged_sccs.erase(it);
                } else {
                    ++it;
                }
            }
        }
    }

    assert(next_pair.first != -1);
    assert(next_pair.second != -1);
    --remaining_merges;
    return next_pair;
}

string MergeSCCs::name() const {
    return "sccs";
}

static shared_ptr<MergeStrategy>_parse(OptionParser &parser) {
    vector<string> order_of_sccs;
    order_of_sccs.push_back("topological");
    order_of_sccs.push_back("reverse_topological");
    order_of_sccs.push_back("decreasing");
    order_of_sccs.push_back("increasing");
    parser.add_enum_option("order_of_sccs",
                           order_of_sccs,
                           "choose an ordering of the sccs: linear (specify "
                           "variable_order) or dfp (specify dfp order options).",
                           "topological");
    vector<string> internal_merge_order;
    internal_merge_order.push_back("linear");
    internal_merge_order.push_back("dfp");
    internal_merge_order.push_back("linear_manual");
    parser.add_enum_option("internal_merge_order",
                           internal_merge_order,
                           "choose an internal merge order: linear (specify "
                           "variable_order) or dfp (specify dfp order options).",
                           "dfp");
    vector<string> merged_sccs_merge_order;
    merged_sccs_merge_order.push_back("linear");
    merged_sccs_merge_order.push_back("dfp");
    parser.add_enum_option("merged_sccs_merge_order",
                           merged_sccs_merge_order,
                           "choose an ordering of the sccs",
                           "dfp");
    // linear merge strategy option
    vector<string> variable_order;
    variable_order.push_back("CG_GOAL_LEVEL");
    variable_order.push_back("CG_GOAL_RANDOM");
    variable_order.push_back("GOAL_CG_LEVEL");
    variable_order.push_back("RANDOM");
    variable_order.push_back("LEVEL");
    variable_order.push_back("REVERSE_LEVEL");
    parser.add_enum_option("variable_order",
                           variable_order,
                           "option useful if merge_order = linear. "
                           "see VariableOrderFinder",
                           "reverse_level");
    // dfp merge strategy options
    vector<string> atomic_ts_order;
    atomic_ts_order.push_back("REGULAR");
    atomic_ts_order.push_back("INVERSE");
    atomic_ts_order.push_back("RANDOM");
    parser.add_enum_option("atomic_ts_order",
                           atomic_ts_order,
                           "order of atomic transition systems",
                           "REGULAR");
    vector<string> product_ts_order;
    product_ts_order.push_back("OLD_TO_NEW");
    product_ts_order.push_back("NEW_TO_OLD");
    product_ts_order.push_back("RANDOM");
    parser.add_enum_option("product_ts_order",
                           product_ts_order,
                           "order of product transition systems",
                           "NEW_TO_OLD");
    parser.add_option<bool>("atomic_before_product",
                            "atomic ts before product ts",
                            "false");
    parser.add_option<bool>("randomized_order",
                            "globally randomized order",
                            "false");
    Options options = parser.parse();
    if (parser.dry_run())
        return 0;
    else
        return make_shared<MergeSCCs>(options);
}

static PluginShared<MergeStrategy> _plugin("merge_sccs", _parse);
