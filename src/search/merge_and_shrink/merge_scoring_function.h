#ifndef MERGE_AND_SHRINK_MERGE_SCORING_FUNCTION_H
#define MERGE_AND_SHRINK_MERGE_SCORING_FUNCTION_H

#include <memory>
#include <vector>

class AbstractTask;

namespace merge_and_shrink {
class FactoredTransitionSystem;
class MergeScoringFunction {
protected:
    bool initialized;
public: // for statistics
    virtual std::string name() const = 0;
protected:
    virtual void dump_specific_options() const {}
public:
    MergeScoringFunction();
    virtual std::vector<int> compute_scores(
        FactoredTransitionSystem &fts,
        const std::vector<std::pair<int, int>> &merge_candidates) = 0;
    virtual void initialize(std::shared_ptr<AbstractTask>) {}
    void dump_options() const;
};
}

#endif
