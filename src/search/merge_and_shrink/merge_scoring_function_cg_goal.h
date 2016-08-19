#ifndef MERGE_AND_SHRINK_MERGE_SCORING_FUNCTION_CG_GOAL_H
#define MERGE_AND_SHRINK_MERGE_SCORING_FUNCTION_CG_GOAL_H

#include "merge_scoring_function.h"

namespace options {
class Options;
}

namespace merge_and_shrink {
class MergeScoringFunctionCgGoal : public MergeScoringFunction {
    bool cg_before_goal;
    std::shared_ptr<AbstractTask> task;
    std::vector<bool> is_goal_variable;
protected:
    virtual std::string name() const override;
public:
    explicit MergeScoringFunctionCgGoal(const options::Options &options);
    virtual ~MergeScoringFunctionCgGoal() override = default;
    virtual std::vector<double> compute_scores(
        FactoredTransitionSystem &fts,
        const std::vector<std::pair<int, int>> &merge_candidates) override;
    virtual void initialize(const std::shared_ptr<AbstractTask> &task) override;
};
}

#endif
