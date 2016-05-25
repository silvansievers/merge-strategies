#ifndef MERGE_AND_SHRINK_MERGE_STRATEGY_FACTORY_DYNAMIC_WEIGHTED_H
#define MERGE_AND_SHRINK_MERGE_STRATEGY_FACTORY_DYNAMIC_WEIGHTED_H

#include "merge_strategy_factory.h"

#include "../options/options.h"

namespace merge_and_shrink {
class Features;

class MergeStrategyFactoryWeighted : public MergeStrategyFactory {
    const options::Options options;
    std::unique_ptr<Features> features;

    std::vector<int> compute_ts_order(
        std::shared_ptr<AbstractTask> task,
        const options::Options &options);

    virtual void dump_strategy_specific_options() const override;
public:
    explicit MergeStrategyFactoryWeighted(const options::Options &opts);
    virtual ~MergeStrategyFactoryWeighted() override = default;
    std::unique_ptr<MergeStrategy> compute_merge_strategy(
        std::shared_ptr<AbstractTask> task,
        FactoredTransitionSystem &fts) override;
    virtual std::string name() const override;
};
}

#endif
