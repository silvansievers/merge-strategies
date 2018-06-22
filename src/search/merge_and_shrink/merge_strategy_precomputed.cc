#include "merge_strategy_precomputed.h"

#include "factored_transition_system.h"
#include "merge_tree.h"

#include "../utils/system.h"

#include <cassert>
#include <iostream>

using namespace std;

namespace merge_and_shrink {
MergeStrategyPrecomputed::MergeStrategyPrecomputed(
    const FactoredTransitionSystem &fts, unique_ptr<MergeTree> merge_tree)
    : MergeStrategy(fts), merge_tree(move(merge_tree)) {
}

pair<int, int> MergeStrategyPrecomputed::get_next(
    const vector<int> &allowed_indices) {
    if (!allowed_indices.empty()) {
        cerr << "Precomputed merge strategies are not compatible with being "
                "computed for a subset of indices" << endl;
        utils::exit_with(utils::ExitCode::UNSUPPORTED);
    }
    assert(!merge_tree->done());
    int next_merge_index = fts.get_size();
    pair<int, int> next_merge = merge_tree->get_next_merge(next_merge_index);
    assert(fts.is_active(next_merge.first));
    assert(fts.is_active(next_merge.second));
    return next_merge;
}
}
