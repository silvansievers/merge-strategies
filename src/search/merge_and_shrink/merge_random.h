#ifndef MERGE_AND_SHRINK_MERGE_RANDOM_H
#define MERGE_AND_SHRINK_MERGE_RANDOM_H

#include "merge_strategy.h"

namespace options {
class Options;
}
namespace utils {
class RandomNumberGenerator;
}

namespace merge_and_shrink {
class MergeRandom : public MergeStrategy {
    const int random_seed;
    std::unique_ptr<utils::RandomNumberGenerator> rng;
protected:
    virtual void dump_strategy_specific_options() const override;
public:
    explicit MergeRandom(const options::Options &options);
    virtual ~MergeRandom() override = default;

    virtual std::pair<int, int> get_next(
        FactoredTransitionSystem &fts) override;
    virtual std::string name() const override;
};
}

#endif
