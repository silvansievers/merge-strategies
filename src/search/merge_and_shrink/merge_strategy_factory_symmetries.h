#ifndef MERGE_AND_SHRINK_MERGE_STRATEGY_FACTORY_SYMMETRIES_H
#define MERGE_AND_SHRINK_MERGE_STRATEGY_FACTORY_SYMMETRIES_H

#include "merge_strategy_factory.h"

#include "../options/options.h"

namespace merge_and_shrink {
class MergeDFP;
class MiasmMergeTree;

class MergeStrategyFactorySymmetries : public MergeStrategyFactory {
    options::Options options;
protected:
    virtual void dump_strategy_specific_options() const override;
public:
    explicit MergeStrategyFactorySymmetries(const options::Options &options);
    virtual ~MergeStrategyFactorySymmetries() override = default;
    virtual std::unique_ptr<MergeStrategy> compute_merge_strategy(
        std::shared_ptr<AbstractTask> task,
        FactoredTransitionSystem &fts) override;
    virtual std::string name() const;
};
}

#endif
