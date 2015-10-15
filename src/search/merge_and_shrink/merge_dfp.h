#ifndef MERGE_AND_SHRINK_MERGE_DFP_H
#define MERGE_AND_SHRINK_MERGE_DFP_H

#include "merge_strategy.h"

class Options;

class MergeDFP : public MergeStrategy {
    bool use_cost;
    // border_atomic_composites is the first index at which a composite
    // transition system can be found in vector of all transition systems as passed
    // as argument to the get_next method.
    int border_atomics_composites;
    int get_corrected_index(int index) const;
    void compute_label_ranks(const TransitionSystem *transition_system,
                             std::vector<int> &label_ranks) const;
    void compute_relevant_labels(const TransitionSystem *ts,
                                 std::vector<bool> &relevant_labels) const;
    std::pair<int, int> get_next_cost(const std::vector<TransitionSystem *> &all_transition_systems);
protected:
    virtual void dump_strategy_specific_options() const override {}
public:
    explicit MergeDFP(const Options &options);
    virtual ~MergeDFP() override = default;
    virtual void initialize(const std::shared_ptr<AbstractTask> task) override;

    virtual std::pair<int, int> get_next(const std::vector<TransitionSystem *> &all_transition_systems) override;
    virtual std::string name() const override;
};

#endif
