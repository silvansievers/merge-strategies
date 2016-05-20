#ifndef MERGE_AND_SHRINK_MERGE_RANDOM_LINEAR_H
#define MERGE_AND_SHRINK_MERGE_RANDOM_LINEAR_H

#include "merge_strategy.h"

#include <memory>
#include <vector>

namespace utils {
class RandomNumberGenerator;
}

namespace merge_and_shrink {
class MergeRandomLinear : public MergeStrategy {
    std::vector<int> randomized_variable_order;
    bool need_first_index;
public:
    MergeRandomLinear(FactoredTransitionSystem &fts,
        std::vector<int> &&randomized_variable_order);
    virtual ~MergeRandomLinear() override = default;
    virtual std::pair<int, int> get_next() override;
};
}

#endif
