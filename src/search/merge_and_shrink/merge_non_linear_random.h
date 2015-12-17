#ifndef MERGE_AND_SHRINK_MERGE_LINEAR_H
#define MERGE_AND_SHRINK_MERGE_LINEAR_H

#include "merge_strategy.h"

class Options;
namespace Utils {
    class RandomNumberGenerator;
}

namespace MergeAndShrink {
class MergeNonLinearRandom : public MergeStrategy {
    const int random_seed;
    std::unique_ptr<Utils::RandomNumberGenerator> rng;
    const int shrink_threshold;
protected:
    virtual void dump_strategy_specific_options() const override;
public:
    explicit MergeNonLinearRandom(const Options &options);
    virtual ~MergeNonLinearRandom() override = default;

    virtual std::pair<int, int> get_next(
        FactoredTransitionSystem &fts) override;
    virtual std::string name() const override;
};
}

#endif
