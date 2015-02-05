#ifndef MERGE_AND_SHRINK_MERGE_SYMMETRIES_LINEAR_H
#define MERGE_AND_SHRINK_MERGE_SYMMETRIES_LINEAR_H

#include "merge_strategy.h"

#include "../option_parser.h"

class MergeSymmetriesLinear : public MergeStrategy {
    // options
    Options options;
    int max_bliss_iterations;
    int bliss_call_time_limit;
    double bliss_remaining_time_budget;

    // statistics
    int iteration_counter;
    int number_of_applied_symmetries;
    bool bliss_limit_reached;
    std::vector<double> bliss_times;
    bool only_applied_dfp;

    std::vector<std::pair<int, int> > merge_order; // TODO: change to from last to first?

    std::vector<int> linear_merge_order;

    void dump_statistics();
protected:
    virtual void dump_strategy_specific_options() const;
public:
    explicit MergeSymmetriesLinear(const Options &options);
    virtual ~MergeSymmetriesLinear() {}

    virtual std::pair<int, int> get_next(const std::vector<TransitionSystem *> &all_abstractions);
    virtual std::string name() const;
};

#endif
