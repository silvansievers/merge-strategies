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
class MergeSelector;
class MergeTree;
class MergeSymmetries : public MergeStrategy {
    options::Options options;
    int num_merges;
    std::shared_ptr<MergeTree> merge_tree;
    std::shared_ptr<MergeSelector> merge_selector;
    int max_bliss_iterations;
    int bliss_call_time_limit;
    double bliss_remaining_time_budget;
    bool tree_is_miasm;

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
