#ifndef MERGE_AND_SHRINK_MERGE_PREDEFINED_H
#define MERGE_AND_SHRINK_MERGE_PREDEFINED_H

#include "merge_strategy.h"

#include <vector>

namespace merge_and_shrink {
class MergePredefined : public MergeStrategy {
    std::vector<std::vector<int>> merge_order;
public:
    MergePredefined(
        FactoredTransitionSystem &fts,
        std::vector<std::vector<int>> merge_order);
    virtual ~MergePredefined() override = default;
    virtual std::pair<int, int> get_next() override;
};
}

#endif
