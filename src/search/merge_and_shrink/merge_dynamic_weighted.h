#ifndef MERGE_AND_SHRINK_MERGE_DYNAMIC_WEIGHTED_H
#define MERGE_AND_SHRINK_MERGE_DYNAMIC_WEIGHTED_H

#include "merge_strategy.h"

#include "../utilities_hash.h"

#include <unordered_map>

class CausalGraph;
class Options;

class MergeDynamicWeighted : public MergeStrategy {
    // Options
    bool debug;
    int w_prefer_causally_connected_vars;
    int w_avoid_additive_vars;
    int w_prefer_small_transitions_states_quotient;
    int w_high_initial_h_value;
    int w_high_average_h_value;
    int w_prefer_ts_large_num_states;
    int w_prefer_ts_large_num_edges;

    // Precomputed stuff
    std::shared_ptr<AbstractTask> task;
    CausalGraph *causal_graph;
    std::vector<int> var_no_to_ts_index;
    std::vector<std::vector<bool> > additive_var_pairs;
    std::unordered_map<std::pair<TransitionSystem *, TransitionSystem *>, double> precomputed_quotients;
    double highest_quotient;
    double lowest_quotient;
    std::unordered_map<std::pair<TransitionSystem *, TransitionSystem *>, double> precomputed_initial_h_improvement;
    double highest_initial_h_improvement;
    double lowest_initial_h_improvement;
    std::unordered_map<std::pair<TransitionSystem *, TransitionSystem *>, double> precomputed_average_h_improvement;
    double highest_average_h_improvement;
    double lowest_average_h_improvement;

    // Statistics
    std::vector<std::pair<int, int> > merge_order;

    // Computation of features
    double compute_feature_causal_connection(
        TransitionSystem *ts1, TransitionSystem *ts2) const;
    double compute_feature_additive_variables(
        TransitionSystem *ts1, TransitionSystem *ts2) const;
    int compute_number_of_product_transitions(
        const TransitionSystem *ts1, const TransitionSystem *ts2) const;
    double compute_feature_transitions_states_quotient(
        TransitionSystem *ts1, TransitionSystem *ts2) const;

    int get_num_transitions(TransitionSystem *ts) const;
    double compute_average_h_value(TransitionSystem *ts) const;
    double normalize_value(double min, double max, double value) const;
    int compute_weighted_sum(
        TransitionSystem *ts1, TransitionSystem *ts2) const;
    void precompute_features(const std::vector<TransitionSystem *> &all_transition_systems);

    virtual void dump_strategy_specific_options() const override;
public:
    MergeDynamicWeighted(const Options opts);
    virtual ~MergeDynamicWeighted();
    virtual void initialize(const std::shared_ptr<AbstractTask> task) override;
    virtual std::pair<int, int> get_next(const std::vector<TransitionSystem *> &all_transition_systems) override;
    virtual std::string name() const override;
};

#endif
