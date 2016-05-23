#ifndef MERGE_AND_SHRINK_MERGE_STRATEGY_H
#define MERGE_AND_SHRINK_MERGE_STRATEGY_H

#include <utility>

namespace merge_and_shrink {
class FactoredTransitionSystem;

class MergeStrategy {
protected:
    FactoredTransitionSystem &fts;
    mutable int iterations_with_tiebreaking;
    mutable int total_tiebreaking_pair_count;
public:
    explicit MergeStrategy(FactoredTransitionSystem &fts);
    virtual ~MergeStrategy() = default;
    virtual std::pair<int, int> get_next() = 0;

    // statistics
    virtual int get_iterations_with_tiebreaking() const {
        return iterations_with_tiebreaking;
    }
    virtual int get_total_tiebreaking_pair_count() const {
        return total_tiebreaking_pair_count;
    }
};
}

#endif
