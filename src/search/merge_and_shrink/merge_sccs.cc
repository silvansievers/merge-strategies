#include "merge_sccs.h"

#include "factored_transition_system.h"
#include "merge_linear.h"
#include "merge_dfp.h"
#include "transition_system.h"

#include "../causal_graph.h"
#include "../option_parser.h"
#include "../plugin.h"
#include "../scc.h"
#include "../variable_order_finder.h"

#include "../utils/logging.h"

#include <algorithm>
#include <cassert>
#include <iostream>

using namespace std;

namespace merge_and_shrink {
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
      merge_dfp(nullptr) {
    options = new Options(options_);
}

MergeSCCs::~MergeSCCs() {
    delete merge_dfp;
    delete options;
}

void MergeSCCs::initialize(const std::shared_ptr<AbstractTask> task) {
    MergeStrategy::initialize(task);

    TaskProxy task_proxy(*task);
    VariablesProxy vars = task_proxy.get_variables();
    int num_vars = vars.size();

    if (internal_merge_order == InternalMergeOrder::DFP || merged_sccs_merge_order == MergedSCCsMergeOrder::DFP) {
        merge_dfp = new MergeDFP(*options);
        merge_dfp->initialize(task);
    }
    if (internal_merge_order == InternalMergeOrder::LINEAR || merged_sccs_merge_order == MergedSCCsMergeOrder::LINEAR) {
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
    for (VariableProxy var : vars) {
        const std::vector<int> &successors =
            task_proxy.get_causal_graph().get_successors(var.get_id());
        cg.push_back(successors);
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
    for (const vector<int> &scc : sccs) {
        cout << scc << endl;
        int scc_size = scc.size();
        if (scc_size == 1) {
            indices_of_merged_sccs.push_back(scc.front());
        } else {
            index += scc_size - 1;
            indices_of_merged_sccs.push_back(index);
            // only store non-singleton sccs for internal merging
            non_singleton_cg_sccs.push_back(vector<int>(scc.begin(), scc.end()));
        }
    }
    if (sccs.size() == 1) {
        cout << "Only one single SCC" << endl;
    }
    if (static_cast<int>(sccs.size()) == num_vars) {
        cout << "Only singleton SCCs" << endl;
        assert(non_singleton_cg_sccs.empty());
    }
//    cout << "indices of merged sccs: " << indices_of_merged_sccs << endl;
}

pair<int, int> MergeSCCs::get_next_linear(
    const FactoredTransitionSystem &fts,
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
                    fts.get_ts(index).get_incorporated_variables();
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
    FactoredTransitionSystem &fts) {
    assert(!done());

    pair<int, int > next_pair = make_pair(-1, -1);
    int most_recent_index = fts.get_size() - 1;
    bool first_merge = false; // needed for linear merging
    if (current_ts_indices.empty()) {
        first_merge = true;
        if (non_singleton_cg_sccs.empty()) {
            assert(indices_of_merged_sccs.size() > 1);
            current_ts_indices = move(indices_of_merged_sccs);
        } else {
            vector<int> &current_scc = non_singleton_cg_sccs.front();
            assert(current_scc.size() > 1);
            current_ts_indices = move(current_scc);
            non_singleton_cg_sccs.erase(non_singleton_cg_sccs.begin());
        }
    } else {
        // Add the newest transition system to the set of current ones
        current_ts_indices.push_back(most_recent_index);
    }

    if (current_ts_indices.size() == 2) {
        next_pair = make_pair(current_ts_indices[0],
                              current_ts_indices[1]);
        current_ts_indices.clear();
    } else {
        if (internal_merge_order == InternalMergeOrder::LINEAR) {
            next_pair = get_next_linear(fts,
                                        current_ts_indices,
                                        most_recent_index,
                                        first_merge);
        } else if (internal_merge_order == InternalMergeOrder::DFP) {
            next_pair = merge_dfp->get_next(fts, current_ts_indices);
        }

        // Remove the two merged indices from the current set of indices
        for (vector<int>::iterator it = current_ts_indices.begin();
             it != current_ts_indices.end();) {
            if (*it == next_pair.first || *it == next_pair.second) {
                it = current_ts_indices.erase(it);
            } else {
                ++it;
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

int MergeSCCs::get_iterations_with_tiebreaking() const {
    if (merge_dfp) {
        return merge_dfp->get_iterations_with_tiebreaking();
    } else {
        return 0;
    }
}

int MergeSCCs::get_total_tiebreaking_pair_count() const {
    if (merge_dfp) {
        return merge_dfp->get_total_tiebreaking_pair_count();
    } else {
        return 0;
    }
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
    MergeLinear::add_options_to_parser(parser);
    // dfp merge strategy options
    MergeDFP::add_options_to_parser(parser);
    Options options = parser.parse();
    if (parser.dry_run())
        return 0;
    else
        return make_shared<MergeSCCs>(options);
}

static PluginShared<MergeStrategy> _plugin("merge_sccs", _parse);
}
