#ifndef MERGE_AND_SHRINK_MERGE_SCORING_FUNCTION_CAUSALLY_CONNECTED_VARIABLE_H
#define MERGE_AND_SHRINK_MERGE_SCORING_FUNCTION_CAUSALLY_CONNECTED_VARIABLE_H

#include "merge_scoring_function.h"

namespace merge_and_shrink {
class MergeScoringFunctionCausallyConnectedVariable : public MergeScoringFunction {
    std::shared_ptr<AbstractTask> task;
protected:
    virtual std::string name() const override;
public:
    MergeScoringFunctionCausallyConnectedVariable() = default;
    virtual ~MergeScoringFunctionCausallyConnectedVariable() override = default;
    virtual std::vector<double> compute_scores(
        FactoredTransitionSystem &fts,
        const std::vector<std::pair<int, int>> &merge_candidates) override;
    virtual void initialize(std::shared_ptr<AbstractTask> task) override;
};
}

#endif
