#ifndef MERGE_AND_SHRINK_MERGE_DYNAMIC_WEIGHTED_H
#define MERGE_AND_SHRINK_MERGE_DYNAMIC_WEIGHTED_H

#include "merge_strategy.h"

class CausalGraph;
class Options;

class MergeDynamicWeighted : public MergeStrategy {
    // Options
    bool debug;
    int w_prefer_causally_connected_vars;
    int w_avoid_additive_vars;
    int w_high_initial_h_value;
    int w_high_average_h_value;
    int w_prefer_ts_large_num_states;
    int w_prefer_ts_large_num_edges;

    // Precomputed stuff
    std::shared_ptr<AbstractTask> task;
    CausalGraph *causal_graph;
    std::vector<int> var_no_to_ts_index;
    std::vector<std::vector<bool> > additive_var_pairs;

    // Statistics
    std::vector<std::pair<int, int> > merge_order;

    double compute_feature_causal_connection(
        TransitionSystem *ts1, TransitionSystem *ts2) const;
    double compute_feature_additive_variables(
        TransitionSystem *ts1, TransitionSystem *ts2) const;
    int get_num_transitions(TransitionSystem *ts) const;
    double get_average_h_value(TransitionSystem *ts) const;
    int compute_number_of_product_transitions(
        const TransitionSystem *ts1, const TransitionSystem *ts2) const;
    int compute_weighted_sum(
        TransitionSystem *ts1, TransitionSystem *ts2) const;

    virtual void dump_strategy_specific_options() const override;
public:
    MergeDynamicWeighted(const Options opts);
    virtual ~MergeDynamicWeighted();
    virtual void initialize(const std::shared_ptr<AbstractTask> task) override;
    virtual std::pair<int, int> get_next(const std::vector<TransitionSystem *> &all_transition_systems) override;
    virtual std::string name() const override;
};

#endif
