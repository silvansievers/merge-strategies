#ifndef MERGE_AND_SHRINK_MERGE_SYMMETRIES_H
#define MERGE_AND_SHRINK_MERGE_SYMMETRIES_H

#include "merge_strategy.h"

#include <memory>
#include <vector>

namespace options {
class Options;
}

namespace merge_and_shrink {
class MergeSelector;
class MergeTree;
class SymmetryGroup;
class MergeSymmetries : public MergeStrategy {
    const int num_merges;

    enum SymmetriesForMerging {
        NO_MERGING,
        SMALLEST,
        LARGEST
    };
    const SymmetriesForMerging symmetries_for_merging;

    enum InternalMerging {
        LINEAR,
        NON_LINEAR
    };
    const InternalMerging internal_merging;

    const int max_bliss_iterations;
    const int bliss_call_time_limit;
    const bool tree_is_miasm;

    std::unique_ptr<SymmetryGroup> symmetry_group;
    double bliss_remaining_time_budget;
    std::shared_ptr<MergeTree> merge_tree;
    std::shared_ptr<MergeSelector> merge_selector;

    // statistics
    int iteration_counter;
    bool bliss_limit_reached;
    std::vector<double> bliss_times;
    bool pure_fallback_strategy;

    // current merge_order
    std::vector<std::pair<int, int>> merge_order;

    void determine_merge_order();
    void dump_statistics();
public:
    MergeSymmetries(
        const FactoredTransitionSystem &fts,
        const options::Options &options,
        int num_merges,
        std::unique_ptr<MergeTree> merge_tree,
        std::shared_ptr<MergeSelector> merge_selector,
        bool tree_is_miasm);
    virtual ~MergeSymmetries() override = default;
    virtual std::pair<int, int> get_next() override;
};
}

#endif
