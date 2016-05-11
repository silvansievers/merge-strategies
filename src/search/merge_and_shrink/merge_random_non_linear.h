#ifndef MERGE_AND_SHRINK_MERGE_RANDOM_NON_LINEAR_H
#define MERGE_AND_SHRINK_MERGE_RANDOM_NON_LINEAR_H

#include "merge_strategy.h"

namespace options {
class Options;
}
namespace utils {
    class RandomNumberGenerator;
}

namespace merge_and_shrink {
class MergeRandomNonLinear : public MergeStrategy {
    int random_seed;
    std::shared_ptr<utils::RandomNumberGenerator> rng;
    const int shrink_threshold;
protected:
    virtual void dump_strategy_specific_options() const override;
public:
    explicit MergeRandomNonLinear(const options::Options &options);
    virtual ~MergeRandomNonLinear() override = default;

    virtual std::pair<int, int> get_next(
        FactoredTransitionSystem &fts) override;
    virtual std::string name() const override;
};
}

#endif
