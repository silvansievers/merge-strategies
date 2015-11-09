#ifndef MERGE_AND_SHRINK_MERGE_LINEAR_H
#define MERGE_AND_SHRINK_MERGE_LINEAR_H

#include "merge_strategy.h"

class Options;
class RandomNumberGenerator;

class MergePredefined : public MergeStrategy {
    std::vector<std::vector<int>> merge_order;
protected:
    virtual void dump_strategy_specific_options() const override;
public:
    explicit MergePredefined(const Options &options);
    virtual ~MergePredefined() override = default;
    virtual void initialize(const std::shared_ptr<AbstractTask> task);

    virtual std::pair<int, int> get_next(
        const std::vector<TransitionSystem *> &all_transition_systems) override;
    virtual std::string name() const override;
};

#endif
