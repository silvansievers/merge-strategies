#ifndef MERGE_AND_SHRINK_MERGE_SYMMETRIES_H
#define MERGE_AND_SHRINK_MERGE_SYMMETRIES_H

#include "merge_strategy.h"

namespace options {
class Options;
}

namespace merge_and_shrink {
class MergeDFP;
class MiasmMergeTree;

class MergeSymmetries : public MergeStrategy {
    options::Options *options;
    int max_bliss_iterations;
    int bliss_call_time_limit;
    double bliss_remaining_time_budget;

    enum FallbackStrategy {
        LINEAR,
        DFP,
        MIASM
    };
    FallbackStrategy fallback_strategy;

    // statistics
    int iteration_counter;
    int number_of_applied_symmetries;
    bool bliss_limit_reached;
    std::vector<double> bliss_times;
    bool pure_fallback_strategy;

    std::vector<std::pair<int, int> > merge_order; // TODO: change to from last to first?

    std::vector<int> linear_merge_order;
    MergeDFP *merge_dfp;
    MiasmMergeTree *miasm_merge_tree;

    void dump_statistics();
    std::pair<int, int> get_next_miasm(FactoredTransitionSystem &fts);
    void update_miasm_merge_tree(
        FactoredTransitionSystem &fts,
        const std::pair<int, int> &next_merge);
protected:
    virtual void dump_strategy_specific_options() const override;
public:
    explicit MergeSymmetries(const options::Options &options);
    virtual ~MergeSymmetries();
    virtual void initialize(const std::shared_ptr<AbstractTask> task) override;

    virtual std::pair<int, int> get_next(FactoredTransitionSystem &fts) override;
    virtual std::string name() const;
};
}

#endif
