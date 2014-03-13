#ifndef MERGE_AND_SHRINK_SYMMETRIES_MERGE_STRATEGY_H
#define MERGE_AND_SHRINK_SYMMETRIES_MERGE_STRATEGY_H

#include "non_linear_merge_strategy.h"

class SymmetriesMergeStrategy : public NonLinearMergeStrategy {
    virtual void dump_strategy_specific_options() const;
public:
    explicit SymmetriesMergeStrategy(const Options &opts);
    virtual ~SymmetriesMergeStrategy() {}

    virtual bool done() const;
    virtual std::pair<int, int> get_next(const std::vector<Abstraction *> &all_abstractions);
    virtual std::string name() const;
};

#endif
