#ifndef MERGE_AND_SHRINK_MERGE_SYMMETRIES_H
#define MERGE_AND_SHRINK_MERGE_SYMMETRIES_H

#include "merge_dfp.h"

#include "../option_parser.h"

class MergeSymmetries : public MergeDFP {
    // options
    Options options;
    int max_bliss_iterations;
    int bliss_call_time_limit;
    double bliss_remaining_time_budget;

    enum FallbackStrategy {
        LINEAR,
        DFP
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

    void dump_statistics();
protected:
    virtual void dump_strategy_specific_options() const;
public:
    explicit MergeSymmetries(const Options &options);
    virtual ~MergeSymmetries() {}
    virtual void initialize(const std::shared_ptr<AbstractTask> task);

    virtual std::pair<int, int> get_next(std::shared_ptr<FactoredTransitionSystem> fts);
    virtual std::string name() const;
};

#endif
