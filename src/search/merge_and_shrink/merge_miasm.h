#ifndef MERGE_AND_SHRINK_MERGE_MIASM_H
#define MERGE_AND_SHRINK_MERGE_MIASM_H

#include "merge_strategy.h"

#include <vector>

namespace merge_and_shrink {
/**
 * @brief The MIASM merging strategy
 * \nosubgrouping
 */
class MergeMiasm : public MergeStrategy {
public:
    MergeMiasm(
        FactoredTransitionSystem &fts,
        std::vector<std::pair<int, int>> merge_order);
    virtual ~MergeMiasm() override = default;
    /** @brief The list of pairs of indices of intermediate abstractions
     * in the order they are merged in */
    std::vector<std::pair<int, int>> merge_order;
public:
    virtual std::pair<int, int> get_next() override;
};
}

#endif
