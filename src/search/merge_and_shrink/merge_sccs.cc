#include "merge_sccs.h"

#include "factored_transition_system.h"
#include "merge_selector_score_based_filtering.h"
#include "transition_system.h"

#include <algorithm>
#include <cassert>
#include <iostream>

using namespace std;

namespace merge_and_shrink {
MergeSCCs::MergeSCCs(FactoredTransitionSystem &fts,
    InternalMergeOrder internal_merge_order,
    std::vector<int> &&linear_variable_order,
    std::shared_ptr<MergeSelectorScoreBasedFiltering> merge_dfp,
    std::vector<std::vector<int>> &&non_singleton_cg_sccs,
    std::vector<int> &&indices_of_merged_sccs)
    : MergeStrategy(fts),
      internal_merge_order(internal_merge_order),
      linear_variable_order(move(linear_variable_order)),
      dfp_selector(merge_dfp),
      non_singleton_cg_sccs(move(non_singleton_cg_sccs)),
      indices_of_merged_sccs(move(indices_of_merged_sccs)) {
}

pair<int, int> MergeSCCs::get_next_linear(
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

pair<int, int> MergeSCCs::get_next() {
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
            next_pair = get_next_linear(current_ts_indices,
                                        most_recent_index,
                                        first_merge);
        } else if (internal_merge_order == InternalMergeOrder::DFP) {
            next_pair = dfp_selector->select_merge(fts, current_ts_indices);
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
    return next_pair;
}

pair<int, int> MergeSCCs::get_dfp_tiebreaking_statistics() const {
    if (dfp_selector) {
        return dfp_selector->get_dfp_tiebreaking_statistics();
    } else {
        return make_pair(0, 0);
    }
}
}
