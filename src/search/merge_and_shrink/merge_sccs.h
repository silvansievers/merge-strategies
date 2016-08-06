#ifndef MERGE_AND_SHRINK_MERGE_SCCS_H
#define MERGE_AND_SHRINK_MERGE_SCCS_H

#include "merge_strategy.h"

#include <memory>
#include <vector>

class AbstractTask;

namespace merge_and_shrink {
class MergeSelector;
class MergeTreeFactory;
class MergeTree;
class MergeSCCs : public MergeStrategy {
    std::shared_ptr<AbstractTask> task;
    std::shared_ptr<MergeTreeFactory> merge_tree_factory;
    std::shared_ptr<MergeSelector> merge_selector;
    std::vector<std::vector<int>> non_singleton_cg_sccs;
    std::vector<int> indices_of_merged_sccs;

    // Active "merge strategies" while merging a set of indices
    std::unique_ptr<MergeTree> current_merge_tree;
    std::vector<int> current_ts_indices;
public:
    MergeSCCs(
        FactoredTransitionSystem &fts,
        std::shared_ptr<AbstractTask> task,
        std::shared_ptr<MergeTreeFactory> merge_tree_factory,
        std::shared_ptr<MergeSelector> merge_selector,
        std::vector<std::vector<int>> non_singleton_cg_sccs,
        std::vector<int> indices_of_merged_sccs);
    virtual ~MergeSCCs() override;
    virtual std::pair<int, int> get_next() override;
    virtual std::pair<int, int> get_tiebreaking_statistics() const override;
};
}

#endif
