#ifndef MERGE_AND_SHRINK_MERGE_SCORING_FUNCTION_CG_GOAL_H
#define MERGE_AND_SHRINK_MERGE_SCORING_FUNCTION_CG_GOAL_H

#include "merge_scoring_function.h"

#include "../task_proxy.h"

namespace plugins {
class Options;
}

namespace merge_and_shrink {
class MergeScoringFunctionCgGoal : public MergeScoringFunction {
    bool cg_before_goal;
    TaskProxy task_proxy;
    std::vector<bool> is_goal_variable;
protected:
    virtual std::string name() const override;
public:
    explicit MergeScoringFunctionCgGoal(const plugins::Options &options);
    virtual ~MergeScoringFunctionCgGoal() override = default;
    virtual std::vector<double> compute_scores(
        const FactoredTransitionSystem &fts,
        const std::vector<std::pair<int, int>> &merge_candidates) override;
    virtual void initialize(const TaskProxy &task_proxy) override;

    virtual bool requires_init_distances() const override {
        return false;
    }

    virtual bool requires_goal_distances() const override {
        return false;
    }
};
}

#endif
