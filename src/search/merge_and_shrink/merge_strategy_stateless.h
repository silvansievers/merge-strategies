#ifndef MERGE_AND_SHRINK_MERGE_STRATEGY_STATELESS_H
#define MERGE_AND_SHRINK_MERGE_STRATEGY_STATELESS_H

#include "merge_strategy.h"

#include <memory>

namespace merge_and_shrink {
class MergeSelector;
class MergeStrategyStateless : public MergeStrategy {
    const std::shared_ptr<MergeSelector> merge_selector;
public:
    MergeStrategyStateless(
        const FactoredTransitionSystem &fts,
        const std::shared_ptr<MergeSelector> &merge_selector);
    virtual ~MergeStrategyStateless() override = default;
    virtual std::pair<int, int> get_next(
        const std::vector<int> &allowed_indices = std::vector<int>()) override;
    virtual std::pair<int, int> get_tiebreaking_statistics() const override;
};
}

#endif
