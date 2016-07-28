#ifndef MERGE_AND_SHRINK_MERGE_SCORING_FUNCTION_MIASM_H
#define MERGE_AND_SHRINK_MERGE_SCORING_FUNCTION_MIASM_H

#include "merge_scoring_function.h"

namespace options {
class Options;
}

namespace merge_and_shrink {
class TransitionSystem;
class MergeScoringFunctionMIASM : public MergeScoringFunction {
    int max_states; // limit to compute shrink sizes
    std::pair<int, int> compute_shrink_sizes(int size1, int size2) const;
protected:
    virtual std::string name() const override;
public:
    explicit MergeScoringFunctionMIASM(const options::Options &options);
    virtual ~MergeScoringFunctionMIASM() override = default;
    virtual std::vector<int> compute_scores(
        FactoredTransitionSystem &fts,
        const std::vector<std::pair<int, int>> &merge_candidates) override;
};
}

#endif
