#ifndef MERGE_AND_SHRINK_MERGE_SYMMETRIES_H
#define MERGE_AND_SHRINK_MERGE_SYMMETRIES_H

#include "merge_strategy.h"

#include "../utils/logging.h"

#include <memory>
#include <vector>

class TaskProxy;

namespace plugins {
class Options;
}

namespace merge_and_shrink {
class MergeSelector;
class MergeTree;
class MergeTreeFactory;
class SymmetryGroup;

enum SymmetriesForMerging {
    NO_MERGING,
    SMALLEST,
    LARGEST
};

enum InternalMerging {
    LINEAR,
    NON_LINEAR,
    SECONDARY
};

class MergeSymmetries : public MergeStrategy {
    const int num_merges;
    const SymmetriesForMerging symmetries_for_merging;
    const InternalMerging internal_merging;
    const int max_bliss_iterations;
    const int bliss_call_time_limit;
    utils::LogProxy log;

    std::unique_ptr<SymmetryGroup> symmetry_group;
    double bliss_remaining_time_budget;
    const TaskProxy &task_proxy;
    std::shared_ptr<MergeTreeFactory> fallback_merge_tree_factory;
    std::shared_ptr<MergeTree> fallback_merge_tree;
    std::shared_ptr<MergeSelector> fallback_merge_selector;
    bool tree_is_miasm;

    // statistics
    int iteration_counter;
    bool bliss_limit_reached;
    std::vector<double> bliss_times;
    bool pure_fallback_strategy;
    bool merging_for_symmetries;
    bool start_merging_for_symmetries;
    bool done_merging_for_symmetries;

    // current merge_order
    std::vector<std::pair<int, int>> merge_order;
    std::unique_ptr<MergeTree> current_merge_tree;
    std::vector<int> current_ts_indices;

    void determine_merge_order();
    void dump_statistics();
public:
    MergeSymmetries(
        const FactoredTransitionSystem &fts,
        const plugins::Options &options,
        const TaskProxy &task_proxy,
        std::shared_ptr<MergeTreeFactory> merge_tree_factory,
        std::shared_ptr<MergeSelector> merge_selector);
    virtual ~MergeSymmetries() override = default;
    virtual std::pair<int, int> get_next(
        const std::vector<int> &allowed_indices = std::vector<int>()) override;
    virtual bool started_merging_for_symmetries() const override {
        return start_merging_for_symmetries;
    }
    virtual bool ended_merging_for_symmetries() const override {
        return done_merging_for_symmetries;
    }
};
}

#endif
