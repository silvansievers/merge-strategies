#ifndef MERGE_AND_SHRINK_MERGE_DYNAMIC_WEIGHTED_H
#define MERGE_AND_SHRINK_MERGE_DYNAMIC_WEIGHTED_H

#include "merge_strategy.h"

#include "../utilities_hash.h"

#include <unordered_map>

class CausalGraph;
class Options;

const int NUM_FEATURES = 7;

//enum FEATURE_TYPE {
//    CAUSAL_CONNECTIONS,
//    NON_ADDITIVITY,
//    TRANS_STATES_QUOT,
//    INIT_H_IMPR,
//    AVG_H_IMPR,
//    INIT_H_SUM,
//    AVG_H_SUM
//};

class AbstractFeature {
    bool merge_required;
public:
    explicit AbstractFeature(bool requires_merge);
    virtual ~AbstractFeature() {}
    virtual double compute_value(const TransitionSystem *ts1,
                                 const TransitionSystem *ts2,
                                 const TransitionSystem *merge) = 0;
    bool requires_merge() const {
        return merge_required;
    }
    virtual void dump() const = 0;
};

class CausalConnectionFeature : public AbstractFeature {
    const std::shared_ptr<AbstractTask> task;
    const CausalGraph &causal_graph;
public:
    CausalConnectionFeature(bool merge_required,
                            const std::shared_ptr<AbstractTask> task);
    virtual double compute_value(const TransitionSystem *ts1,
                                 const TransitionSystem *ts2,
                                 const TransitionSystem *merge) override;
    virtual void dump() const;
};

class NonAdditivityFeature : public AbstractFeature {
    std::vector<std::vector<bool> > additive_var_pairs;
public:
    NonAdditivityFeature(bool merge_required,
                         const std::shared_ptr<AbstractTask> task);
    virtual double compute_value(const TransitionSystem *ts1,
                                 const TransitionSystem *ts2,
                                 const TransitionSystem *merge) override;
    virtual void dump() const;
};

class Features {
    std::vector<int> weights;
    const std::shared_ptr<AbstractTask> task;
    bool debug;
    std::vector<AbstractFeature *> features;
    std::vector<double> min_values;
    std::vector<double> max_values;
    std::unordered_map<std::pair<TransitionSystem *, TransitionSystem *>,
                       std::vector<double> > unnormalized_values;
    void update_min_max(int feature_no, double value);
    double normalize_value(int feature_no, double value) const;
public:
    Features(std::vector<int> &&weights,
             const std::shared_ptr<AbstractTask> task,
             bool debug);
    ~Features();
    void precompute_unnormalized_values(TransitionSystem *ts1,
                                        TransitionSystem *ts2);
    double compute_weighted_normalized_sum(
        TransitionSystem *ts1, TransitionSystem *ts2) const;
    void clear();
};

class MergeDynamicWeighted : public MergeStrategy {
    // Options
    bool debug;
    int w_causally_connected_vars;
    int w_nonadditive_vars;
//    int w_small_transitions_states_quotient;
//    int w_high_initial_h_value_improvement;
//    int w_high_average_h_value_improvement;
//    int w_high_initial_h_value_sum;
//    int w_high_average_h_value_sum;


    // Precomputed stuff
    std::vector<int> var_no_to_ts_index;
    Features *features;

//    std::unordered_map<std::pair<TransitionSystem *, TransitionSystem *>, double> precomputed_quotients;
//    double highest_quotient;
//    double lowest_quotient;
//    std::unordered_map<std::pair<TransitionSystem *, TransitionSystem *>, double> precomputed_initial_h_improvement;
//    double highest_initial_h_improvement;
//    double lowest_initial_h_improvement;
//    std::unordered_map<std::pair<TransitionSystem *, TransitionSystem *>, double> precomputed_average_h_improvement;
//    double highest_average_h_improvement;
//    double lowest_average_h_improvement;
//    int highest_initial_h_sum;
//    int lowest_initial_h_sum;
//    std::unordered_map<std::pair<TransitionSystem *, TransitionSystem *>, double> precomputed_average_h_sum;
//    double highest_average_h_sum;
//    double lowest_average_h_sum;

    // Statistics
    std::vector<std::pair<int, int> > merge_order;

    // Computation of features
    int compute_number_of_product_transitions(
        const TransitionSystem *ts1, const TransitionSystem *ts2) const;
    double compute_feature_transitions_states_quotient(
        TransitionSystem *ts1, TransitionSystem *ts2) const;
    int get_num_transitions(TransitionSystem *ts) const;
    double compute_average_h_value(TransitionSystem *ts) const;
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
