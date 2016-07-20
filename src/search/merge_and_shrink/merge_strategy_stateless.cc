#include "merge_strategy_stateless.h"

#include "merge_selector.h"
#include "merge_selector_score_based_filtering.h"

using namespace std;

namespace merge_and_shrink {
MergeStrategyStateless::MergeStrategyStateless(
    FactoredTransitionSystem &fts,
    std::shared_ptr<MergeSelector> merge_selector)
    : MergeStrategy(fts),
      merge_selector(merge_selector) {
}

pair<int, int> MergeStrategyStateless::get_next() {
    return merge_selector->select_merge(fts);
}

pair<int, int> MergeStrategyStateless::get_dfp_tiebreaking_statistics() const {
    shared_ptr<MergeSelectorScoreBasedFiltering> cast =
        dynamic_pointer_cast<MergeSelectorScoreBasedFiltering>(merge_selector);
    if (cast) {
        return cast->get_dfp_tiebreaking_statistics();
    } else {
        return make_pair(0, 0);
    }
}
}