#ifndef MERGE_AND_SHRINK_MERGE_SYMMETRIES_H
#define MERGE_AND_SHRINK_MERGE_SYMMETRIES_H

#include "merge_strategy.h"

#include "../options/options.h"

#include <memory>
#include <vector>

namespace options {
class Options;
}

namespace merge_and_shrink {
class MergeSelectorScoreBasedFiltering;
class MergeTree;

enum FallbackStrategy {
    LINEAR,
    DFP,
    MIASM
};

class MergeSymmetries : public MergeStrategy {
    options::Options options;
    int num_merges;
    std::vector<int> linear_merge_order;
    std::shared_ptr<MergeSelectorScoreBasedFiltering> dfp_selector;
    std::unique_ptr<MergeTree> miasm_merge_tree;
    int max_bliss_iterations;
    int bliss_call_time_limit;
    double bliss_remaining_time_budget;
    FallbackStrategy fallback_strategy;

    // statistics
    int iteration_counter;
    int number_of_applied_symmetries;
    bool bliss_limit_reached;
    std::vector<double> bliss_times;
    bool pure_fallback_strategy;

    // current merge_order
    std::vector<std::pair<int, int>> merge_order;

    void dump_statistics();
public:
    MergeSymmetries(
        FactoredTransitionSystem &fts,
        const options::Options &options,
        int num_merges,
        std::vector<int> linear_merge_order,
        std::shared_ptr<MergeSelectorScoreBasedFiltering> dfp_selector,
        std::unique_ptr<MergeTree> miasm_merge_tree);
    virtual ~MergeSymmetries() override = default;
    virtual std::pair<int, int> get_next() override;
};
}

#endif
