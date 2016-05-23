#ifndef MERGE_AND_SHRINK_MERGE_STRATEGY_FACTORY_RANDOM_NON_LINEAR_H
#define MERGE_AND_SHRINK_MERGE_STRATEGY_FACTORY_RANDOM_NON_LINEAR_H

#include "merge_strategy_factory.h"

namespace options {
class Options;
}
namespace utils {
class RandomNumberGenerator;
}

namespace merge_and_shrink {
class MergeStrategyFactoryRandomNonLinear : public MergeStrategyFactory {
    int shrink_threshold;
    int random_seed;
    std::shared_ptr<utils::RandomNumberGenerator> rng;
protected:
    virtual void dump_strategy_specific_options() const override;
public:
    explicit MergeStrategyFactoryRandomNonLinear(const options::Options &options);
    virtual ~MergeStrategyFactoryRandomNonLinear() override = default;
    virtual std::unique_ptr<MergeStrategy> compute_merge_strategy(
        const std::shared_ptr<AbstractTask> &task,
        FactoredTransitionSystem &fts) override;
    virtual std::string name() const override;
};
}

#endif
