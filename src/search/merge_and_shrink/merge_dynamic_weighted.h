#ifndef MERGE_AND_SHRINK_MERGE_DYNAMIC_WEIGHTED_H
#define MERGE_AND_SHRINK_MERGE_DYNAMIC_WEIGHTED_H

#include "merge_strategy.h"

class Options;

class MergeDynamicWeighted : public MergeStrategy {
    int w_prefer_causally_connected_vars;
    int w_avoid_additive_vars;
    int w_high_initial_h_value;
    int w_high_average_h_value;
    int w_prefer_ts_large_num_states;
    int w_prefer_ts_large_num_edges;
    std::shared_ptr<AbstractTask> task;

    std::vector<int> var_no_to_ts_index;
    std::vector<std::vector<bool> > additive_var_pairs;

    virtual void dump_strategy_specific_options() const override;
public:
    MergeDynamicWeighted(const Options opts);
    virtual ~MergeDynamicWeighted() = default;
    virtual void initialize(const std::shared_ptr<AbstractTask> task) override;
    virtual std::pair<int, int> get_next(const std::vector<TransitionSystem *> &all_transition_systems) override;
    virtual std::string name() const override;
};

#endif
