#ifndef MERGE_AND_SHRINK_MERGE_SELECTOR_SCORE_BASED_WEIGHTED_SUM_H
#define MERGE_AND_SHRINK_MERGE_SELECTOR_SCORE_BASED_WEIGHTED_SUM_H

#include "merge_selector.h"

#include "merge_scoring_function.h"

#include <vector>

namespace options {
class Options;
}

namespace merge_and_shrink {
extern const int MINUSINF;
extern double normalize_value(double min_score, double max_score, double score);

class MergeSelectorScoreBasedWeightedSum : public MergeSelector {
    std::vector<std::shared_ptr<MergeScoringFunction>> merge_scoring_functions;
    std::vector<int> weights;
    bool normalize;
protected:
    virtual std::string name() const override;
    virtual void dump_specific_options() const override;
public:
    explicit MergeSelectorScoreBasedWeightedSum(const options::Options &options);
    virtual ~MergeSelectorScoreBasedWeightedSum() override = default;
    virtual std::pair<int, int> select_merge(
        FactoredTransitionSystem &fts) const override;
    virtual void initialize(std::shared_ptr<AbstractTask> task) override;
};
}

#endif
