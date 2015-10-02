#ifndef MERGE_AND_SHRINK_MERGE_DYNAMIC_WEIGHTED_H
#define MERGE_AND_SHRINK_MERGE_DYNAMIC_WEIGHTED_H

#include "merge_strategy.h"

#include "../utilities_hash.h"

#include <unordered_map>

class CausalGraph;
class Options;

const int NUM_FEATURES = 7;

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
    virtual void dump_precomputed_data() const = 0;
};

class CausalConnectionFeature : public AbstractFeature {
    const std::shared_ptr<AbstractTask> task;
    const CausalGraph &causal_graph;
public:
    explicit CausalConnectionFeature(const std::shared_ptr<AbstractTask> task);
    virtual double compute_value(const TransitionSystem *ts1,
                                 const TransitionSystem *ts2,
                                 const TransitionSystem *merge) override;
    virtual void dump_precomputed_data() const;
};

class NonAdditivityFeature : public AbstractFeature {
    std::vector<std::vector<bool> > additive_var_pairs;
public:
    explicit NonAdditivityFeature(const std::shared_ptr<AbstractTask> task);
    virtual double compute_value(const TransitionSystem *ts1,
                                 const TransitionSystem *ts2,
                                 const TransitionSystem *merge) override;
    virtual void dump_precomputed_data() const;
};

class TransStatesQuotFeature : public AbstractFeature {
public:
    TransStatesQuotFeature();
    virtual double compute_value(const TransitionSystem *ts1,
                                 const TransitionSystem *ts2,
                                 const TransitionSystem *merge) override;
    virtual void dump_precomputed_data() const {}
};

class InitHImprovementFeature : public AbstractFeature {
public:
    InitHImprovementFeature();
    virtual double compute_value(const TransitionSystem *ts1,
                                 const TransitionSystem *ts2,
                                 const TransitionSystem *merge) override;
    virtual void dump_precomputed_data() const {}
};

class AvgHImprovementFeature : public AbstractFeature {
public:
    AvgHImprovementFeature();
    virtual double compute_value(const TransitionSystem *ts1,
                                 const TransitionSystem *ts2,
                                 const TransitionSystem *merge) override;
    virtual void dump_precomputed_data() const {}
};

class InitHSumFeature : public AbstractFeature {
public:
    InitHSumFeature();
    virtual double compute_value(const TransitionSystem *ts1,
                                 const TransitionSystem *ts2,
                                 const TransitionSystem *merge) override;
    virtual void dump_precomputed_data() const {}
};

class AvgHSumFeature : public AbstractFeature {
public:
    AvgHSumFeature();
    virtual double compute_value(const TransitionSystem *ts1,
                                 const TransitionSystem *ts2,
                                 const TransitionSystem *merge) override;
    virtual void dump_precomputed_data() const {}
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
    // TODO: move weights inside features. construct features in constructor,
    // and give it an initialize method
    int w_causally_connected_vars;
    int w_nonadditive_vars;
    int w_small_transitions_states_quotient;
    int w_high_initial_h_value_improvement;
    int w_high_average_h_value_improvement;
    int w_high_initial_h_value_sum;
    int w_high_average_h_value_sum;

    // Precomputed stuff
    std::vector<int> var_no_to_ts_index;
    Features *features;

    // Statistics
    std::vector<std::pair<int, int> > merge_order;

    virtual void dump_strategy_specific_options() const override;
public:
    MergeDynamicWeighted(const Options opts);
    virtual ~MergeDynamicWeighted();
    virtual void initialize(const std::shared_ptr<AbstractTask> task) override;
    virtual std::pair<int, int> get_next(const std::vector<TransitionSystem *> &all_transition_systems) override;
    virtual std::string name() const override;
};

#endif
