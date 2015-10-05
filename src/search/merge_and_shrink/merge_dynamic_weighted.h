#ifndef MERGE_AND_SHRINK_MERGE_DYNAMIC_WEIGHTED_H
#define MERGE_AND_SHRINK_MERGE_DYNAMIC_WEIGHTED_H

#include "merge_strategy.h"

#include "../utilities_hash.h"

#include <unordered_map>

class CausalGraph;
class Options;
class TaskProxy;

const int NUM_FEATURES = 7;

class Feature {
    const int id;
    const std::string name;
    const bool merge_required;
    const int weight;
    virtual double compute_value(const TransitionSystem *ts1,
                                 const TransitionSystem *ts2,
                                 const TransitionSystem *merge) = 0;
public:
    Feature(int id, std::string name, bool requires_merge, int weight);
    virtual ~Feature() {}
    double compute_unnormalized_value(const TransitionSystem *ts1,
                                      const TransitionSystem *ts2,
                                      const TransitionSystem *merge);
    int get_id() const {
        return id;
    }
    std::string get_name() const {
        return name;
    }
    bool requires_merge() const {
        return merge_required;
    }
    int get_weight() const {
        return weight;
    }
    virtual void initialize(const TaskProxy &, bool) {}
    void dump() const;
};

class CausalConnectionFeature : public Feature {
    CausalGraph *causal_graph;
    virtual double compute_value(const TransitionSystem *ts1,
                                 const TransitionSystem *ts2,
                                 const TransitionSystem *merge) override;
public:
    CausalConnectionFeature(int id, int weight);
    virtual ~CausalConnectionFeature();
    virtual void initialize(const TaskProxy &task_proxy, bool dump) override;
};

class NonAdditivityFeature : public Feature {
    std::vector<std::vector<bool> > additive_var_pairs;
    virtual double compute_value(const TransitionSystem *ts1,
                                 const TransitionSystem *ts2,
                                 const TransitionSystem *merge) override;
public:
    NonAdditivityFeature(int id, int weight);
    virtual void initialize(const TaskProxy &task_proxy, bool dump) override;
};

class TransStatesQuotFeature : public Feature {
    virtual double compute_value(const TransitionSystem *ts1,
                                 const TransitionSystem *ts2,
                                 const TransitionSystem *merge) override;
public:
    TransStatesQuotFeature(int id, int weight);
};

class InitHImprovementFeature : public Feature {
    virtual double compute_value(const TransitionSystem *ts1,
                                 const TransitionSystem *ts2,
                                 const TransitionSystem *merge) override;
public:
    InitHImprovementFeature(int id, int weight);
};

class AvgHImprovementFeature : public Feature {
    virtual double compute_value(const TransitionSystem *ts1,
                                 const TransitionSystem *ts2,
                                 const TransitionSystem *merge) override;
public:
    AvgHImprovementFeature(int id, int weight);
};

class InitHSumFeature : public Feature {
    virtual double compute_value(const TransitionSystem *ts1,
                                 const TransitionSystem *ts2,
                                 const TransitionSystem *merge) override;
public:
    InitHSumFeature(int id, int weight);
};

class AvgHSumFeature : public Feature {
    virtual double compute_value(const TransitionSystem *ts1,
                                 const TransitionSystem *ts2,
                                 const TransitionSystem *merge) override;
public:
    AvgHSumFeature(int id, int weight);
};

class Features {
    const bool debug;
    TaskProxy *task_proxy;
    std::vector<Feature *> features;
    std::vector<double> min_values;
    std::vector<double> max_values;
    std::unordered_map<std::pair<TransitionSystem *, TransitionSystem *>,
                       std::vector<double> > unnormalized_values;
    void update_min_max(int feature_no, double value);
    double normalize_value(int feature_no, double value) const;
public:
    explicit Features(const Options opts);
    ~Features();
    void initialize(const std::shared_ptr<AbstractTask> task);
    void precompute_unnormalized_values(TransitionSystem *ts1,
                                        TransitionSystem *ts2);
    double compute_weighted_normalized_sum(
        TransitionSystem *ts1, TransitionSystem *ts2) const;
    void clear();
    void dump_weights() const;
};

class MergeDynamicWeighted : public MergeStrategy {
    Features *features;
    std::vector<int> var_no_to_ts_index;

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
