#ifndef MERGE_AND_SHRINK_MERGE_STRATEGY_FACTORY_PREDEFINED_H
#define MERGE_AND_SHRINK_MERGE_STRATEGY_FACTORY_PREDEFINED_H

#include "merge_strategy_factory.h"

#include <vector>

namespace options {
class Options;
}

namespace merge_and_shrink {
class BinaryTree;

class MergeStrategyFactoryPredefined : public MergeStrategyFactory {
    std::vector<std::vector<int>> merge_order;
    BinaryTree *root;
protected:
    virtual void dump_strategy_specific_options() const override;
public:
    explicit MergeStrategyFactoryPredefined(const options::Options &options);
    virtual ~MergeStrategyFactoryPredefined() override;
    virtual std::unique_ptr<MergeStrategy> compute_merge_strategy(
        const std::shared_ptr<AbstractTask> &task,
        FactoredTransitionSystem &fts) override;
    virtual std::string name() const override;
};
}

#endif
