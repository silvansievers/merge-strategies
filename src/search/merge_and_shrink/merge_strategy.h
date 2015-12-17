#ifndef MERGE_AND_SHRINK_MERGE_STRATEGY_H
#define MERGE_AND_SHRINK_MERGE_STRATEGY_H

#include <memory>
#include <string>
#include <utility>
#include <vector>

class AbstractTask;


namespace MergeAndShrink {
class FactoredTransitionSystem;

class MergeStrategy {
    const int UNINITIALIZED = -1;
protected:
    int remaining_merges;
    int iterations_with_tiebreaking;
    int total_tiebreaking_pair_count;
    bool initialized() const;
    virtual void dump_strategy_specific_options() const = 0;
public:
    MergeStrategy();
    virtual ~MergeStrategy() = default;
    virtual void initialize(const std::shared_ptr<AbstractTask> task);
    bool done() const;
    void dump_options() const;

    // Implementations of get_next have to decrease remaining_merges by one
    virtual std::pair<int, int> get_next(FactoredTransitionSystem &fts) = 0;
    virtual std::string name() const = 0;

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
