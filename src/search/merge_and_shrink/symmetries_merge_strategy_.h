#ifndef MERGE_AND_SHRINK_SYMMETRIES_MERGE_STRATEGY_H
#define MERGE_AND_SHRINK_SYMMETRIES_MERGE_STRATEGY_H

#include "non_linear_merge_strategy.h"

class SymmetriesMergeStrategy : public NonLinearMergeStrategy {
public:
    explicit NonLinearMergeStrategy(const Options &opts);
    virtual ~NonLinearMergeStrategy() {}

    virtual bool done() const;
    virtual void get_next(const std::vector<Abstraction *> &all_abstractions, std::pair<int, int> &next_indices);
    virtual std::string name() const;
    virtual void print_summary() const;
};

#endif
