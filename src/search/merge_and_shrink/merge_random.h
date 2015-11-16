#ifndef MERGE_AND_SHRINK_MERGE_LINEAR_H
#define MERGE_AND_SHRINK_MERGE_LINEAR_H

#include "merge_strategy.h"

class Options;
class RandomNumberGenerator;

class MergeRandom : public MergeStrategy {
    const int random_seed;
    std::unique_ptr<RandomNumberGenerator> rng;
protected:
    virtual void dump_strategy_specific_options() const override;
public:
    explicit MergeRandom(const Options &options);
    virtual ~MergeRandom() override = default;

    virtual std::pair<int, int> get_next(
        std::shared_ptr<FactoredTransitionSystem> fts) override;
    virtual std::string name() const override;
};

#endif