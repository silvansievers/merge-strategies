#ifndef MERGE_AND_SHRINK_MERGE_RANDOM_H
#define MERGE_AND_SHRINK_MERGE_RANDOM_H

#include "merge_strategy.h"

#include <memory>

namespace utils {
class RandomNumberGenerator;
}

namespace merge_and_shrink {
class MergeRandom : public MergeStrategy {
    std::shared_ptr<utils::RandomNumberGenerator> rng;
public:
    MergeRandom(FactoredTransitionSystem &fts,
                std::shared_ptr<utils::RandomNumberGenerator> rng);
    virtual ~MergeRandom() override = default;
    virtual std::pair<int, int> get_next() override;
};
}

#endif
