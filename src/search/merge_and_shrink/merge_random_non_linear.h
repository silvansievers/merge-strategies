#ifndef MERGE_AND_SHRINK_MERGE_RANDOM_NON_LINEAR_H
#define MERGE_AND_SHRINK_MERGE_RANDOM_NON_LINEAR_H

#include "merge_strategy.h"

#include <memory>

namespace utils {
class RandomNumberGenerator;
}

namespace merge_and_shrink {
class MergeRandomNonLinear : public MergeStrategy {
    std::shared_ptr<utils::RandomNumberGenerator> rng;
    int shrink_threshold;
public:
    MergeRandomNonLinear(
        FactoredTransitionSystem &fts,
        std::shared_ptr<utils::RandomNumberGenerator> rng,
        int shrink_threshold);
    virtual ~MergeRandomNonLinear() override = default;
    virtual std::pair<int, int> get_next() override;
};
}

#endif
