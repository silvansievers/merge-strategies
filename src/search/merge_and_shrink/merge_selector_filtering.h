#ifndef MERGE_AND_SHRINK_MERGE_SELECTOR_FILTERING_H
#define MERGE_AND_SHRINK_MERGE_SELECTOR_FILTERING_H

#include "merge_selector.h"

#include <memory>
#include <vector>

namespace plugins {
class Options;
}

namespace merge_and_shrink {
struct MergeCandidate;
class MergeScoringFunction;
class MergeSelectorFiltering : public MergeSelector {
    std::vector<std::shared_ptr<MergeScoringFunction>> merge_scoring_functions;
    mutable int iterations_with_tiebreaking;
    mutable int total_tiebreaking_pair_count;
    mutable std::vector<std::vector<std::shared_ptr<MergeCandidate>>> merge_candidates_by_indices;
    mutable int num_candidates;

    std::shared_ptr<MergeCandidate> get_candidate(int index1, int index2) const;
protected:
    virtual std::string name() const override;
    virtual void dump_selector_specific_options(utils::LogProxy &log) const override;
public:
    explicit MergeSelectorFiltering(const plugins::Options &options);
    virtual ~MergeSelectorFiltering() override = default;
    virtual std::pair<int, int> select_merge(
        const FactoredTransitionSystem &fts,
        const std::vector<int> &indices_subset = std::vector<int>()) const override;
    virtual void initialize(const TaskProxy &task_proxy) override;
    virtual bool requires_init_distances() const override;
    virtual bool requires_goal_distances() const override;
    std::pair<int, int> get_tiebreaking_statistics() const {
        return std::make_pair(iterations_with_tiebreaking,
                              total_tiebreaking_pair_count);
    }
};
}

#endif
