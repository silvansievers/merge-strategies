#ifndef MERGE_AND_SHRINK_MERGE_DYNAMIC_WEIGHTED_H
#define MERGE_AND_SHRINK_MERGE_DYNAMIC_WEIGHTED_H

#include "merge_strategy.h"

#include "../utilities_hash.h"

#include <unordered_map>

class Options;
class TaskProxy;

const int NUM_FEATURES = 7;

class Feature {
    const int id;
    const int weight;
    const std::string name;
    const bool merge_required;
    const bool minimize_value;
    virtual double compute_value(
        const std::shared_ptr<FactoredTransitionSystem> fts,
        int ts_index1,
        int ts_index2,
        int merge_index) = 0;
public:
    Feature(int id, int weight, std::string name,
            bool requires_merge, bool minimize_value);
    virtual ~Feature() {}
    virtual void initialize(const TaskProxy &, bool) {}
    virtual void precompute_data(
        const std::shared_ptr<FactoredTransitionSystem>) {}
    double compute_unnormalized_value(
        const std::shared_ptr<FactoredTransitionSystem> fts,
        int ts_index1,
        int ts_index2,
        int merge_index);
    virtual void clear() {}
    int get_id() const {
        return id;
    }
    int get_weight() const {
        return weight;
    }
    std::string get_name() const {
        return name;
    }
    bool requires_merge() const {
        return merge_required;
    }
    bool minimize() const {
        return minimize_value;
    }
    void dump() const;
};

class CausalConnectionFeature : public Feature {
    std::vector<std::vector<int>> var_pair_causal_links;
    virtual double compute_value(
        const std::shared_ptr<FactoredTransitionSystem> fts,
        int ts_index1,
        int ts_index2,
        int merge_index) override;
public:
    CausalConnectionFeature(int id, int weight);
    virtual void initialize(const TaskProxy &task_proxy, bool dump) override;
};

class BoolCausalConnectionFeature : public Feature {
    std::vector<std::vector<bool>> causally_connected_var_pairs;
    virtual double compute_value(
        const std::shared_ptr<FactoredTransitionSystem> fts,
        int ts_index1,
        int ts_index2,
        int merge_index) override;
public:
    BoolCausalConnectionFeature(int id, int weight);
    virtual void initialize(const TaskProxy &task_proxy, bool dump) override;
};

class NonAdditivityFeature : public Feature {
    std::vector<std::vector<bool>> additive_var_pairs;
    virtual double compute_value(
        const std::shared_ptr<FactoredTransitionSystem> fts,
        int ts_index1,
        int ts_index2,
        int merge_index) override;
public:
    NonAdditivityFeature(int id, int weight);
    virtual void initialize(const TaskProxy &task_proxy, bool dump) override;
};

class SmallTransStatesQuotFeature : public Feature {
    virtual double compute_value(
        const std::shared_ptr<FactoredTransitionSystem> fts,
        int ts_index1,
        int ts_index2,
        int merge_index) override;
public:
    SmallTransStatesQuotFeature(int id, int weight);
};

class HighTransStatesQuotFeature : public Feature {
    virtual double compute_value(
        const std::shared_ptr<FactoredTransitionSystem> fts,
        int ts_index1,
        int ts_index2,
        int merge_index) override;
public:
    HighTransStatesQuotFeature(int id, int weight);
};

class InitHImprovementFeature : public Feature {
    virtual double compute_value(
        const std::shared_ptr<FactoredTransitionSystem> fts,
        int ts_index1,
        int ts_index2,
        int merge_index) override;
public:
    InitHImprovementFeature(int id, int weight);
};

class AbsoluteInitHFeature : public Feature {
    virtual double compute_value(
        const std::shared_ptr<FactoredTransitionSystem> fts,
        int ts_index1,
        int ts_index2,
        int merge_index) override;
public:
    AbsoluteInitHFeature(int id, int weight);
};

class AbsoluteMaxFFeature : public Feature {
    virtual double compute_value(
        const std::shared_ptr<FactoredTransitionSystem> fts,
        int ts_index1,
        int ts_index2,
        int merge_index) override;
public:
    AbsoluteMaxFFeature(int id, int weight);
};

class AbsoluteMaxGFeature : public Feature {
    virtual double compute_value(
        const std::shared_ptr<FactoredTransitionSystem> fts,
        int ts_index1,
        int ts_index2,
        int merge_index) override;
public:
    AbsoluteMaxGFeature(int id, int weight);
};

class AbsoluteMaxHFeature : public Feature {
    virtual double compute_value(
        const std::shared_ptr<FactoredTransitionSystem> fts,
        int ts_index1,
        int ts_index2,
        int merge_index) override;
public:
    AbsoluteMaxHFeature(int id, int weight);
};

class AvgHImprovementFeature : public Feature {
    virtual double compute_value(
        const std::shared_ptr<FactoredTransitionSystem> fts,
        int ts_index1,
        int ts_index2,
        int merge_index) override;
public:
    AvgHImprovementFeature(int id, int weight);
};

class InitHSumFeature : public Feature {
    virtual double compute_value(
        const std::shared_ptr<FactoredTransitionSystem> fts,
        int ts_index1,
        int ts_index2,
        int merge_index) override;
public:
    InitHSumFeature(int id, int weight);
};

class AvgHSumFeature : public Feature {
    virtual double compute_value(
        const std::shared_ptr<FactoredTransitionSystem> fts,
        int ts_index1,
        int ts_index2,
        int merge_index) override;
public:
    AvgHSumFeature(int id, int weight);
};

class DFPFeature : public Feature {
    std::vector<std::vector<int>> ts_index_to_label_ranks;
    virtual double compute_value(
        const std::shared_ptr<FactoredTransitionSystem> fts,
        int ts_index1,
        int ts_index2,
        int merge_index) override;
public:
    DFPFeature(int id, int weight);
    virtual void initialize(const TaskProxy &task_proxy, bool dump) override;
//    virtual void precompute_data(const std::shared_ptr<FactoredTransitionSystem> fts) override;
    virtual void clear() override;
};

class GoalRelevanceFeature : public Feature {
    virtual double compute_value(
        const std::shared_ptr<FactoredTransitionSystem> fts,
        int ts_index1,
        int ts_index2,
        int merge_index) override;
public:
    GoalRelevanceFeature(int id, int weight);
};

class NumVariablesFeature : public Feature {
    virtual double compute_value(
        const std::shared_ptr<FactoredTransitionSystem> fts,
        int ts_index1,
        int ts_index2,
        int merge_index) override;
public:
    NumVariablesFeature(int id, int weight);
};

class ShrinkPerfectlyFeature : public Feature {
    virtual double compute_value(
        const std::shared_ptr<FactoredTransitionSystem> fts,
        int ts_index1,
        int ts_index2,
        int merge_index) override;
public:
    ShrinkPerfectlyFeature(int id, int weight);
};

class NumTransitionsFeature : public Feature {
    virtual double compute_value(
        const std::shared_ptr<FactoredTransitionSystem> fts,
        int ts_index1,
        int ts_index2,
        int merge_index) override;
public:
    NumTransitionsFeature(int id, int weight);
};

class LROpportunitiesFeatures : public Feature {
    std::unordered_map<std::pair<int, int>, int> ts_pair_to_combinable_label_count;
    virtual double compute_value(
        const std::shared_ptr<FactoredTransitionSystem> fts,
        int ts_index1,
        int ts_index2,
        int merge_index) override;
public:
    LROpportunitiesFeatures(int id, int weight);
    virtual void clear() override;
    void precompute_data(const std::shared_ptr<FactoredTransitionSystem> fts) override;
};

class MoreLROpportunitiesFeatures : public Feature {
    std::unordered_map<std::pair<int, int>, int> ts_pair_to_combinable_label_count;
    virtual double compute_value(
        const std::shared_ptr<FactoredTransitionSystem> fts,
        int ts_index1,
        int ts_index2,
        int merge_index) override;
public:
    MoreLROpportunitiesFeatures(int id, int weight);
    virtual void clear() override;
    void precompute_data(const std::shared_ptr<FactoredTransitionSystem> fts) override;
};

class MIASMFeature : public Feature {
    virtual double compute_value(
        const std::shared_ptr<FactoredTransitionSystem> fts,
        int ts_index1,
        int ts_index2,
        int merge_index) override;
public:
    MIASMFeature(int id, int weight);
};

class MutexFeature : public Feature {
    virtual double compute_value(
        const std::shared_ptr<FactoredTransitionSystem> fts,
        int ts_index1,
        int ts_index2,
        int merge_index) override;
public:
    MutexFeature(int id, int weight);
};

class Features {
    const bool debug;
    bool merge_required; // if any of the active features require the merge
    std::vector<Feature *> features;
    std::vector<double> min_values; // finite minium values of all features
    std::vector<double> max_values; // finite maximum values of all features
    std::unordered_map<std::pair<int, int>, std::vector<double>> unnormalized_values;
    void update_min_max(int feature_no, double value);
    double normalize_value(int feature_no, double value) const;
public:
    explicit Features(const Options opts);
    ~Features();
    void initialize(const TaskProxy &task_proxy);
    void precompute_data(const std::shared_ptr<FactoredTransitionSystem> fts);
    bool require_merge() const {
        return merge_required;
    }
    void precompute_unnormalized_values(
        const std::shared_ptr<FactoredTransitionSystem> fts,
        int ts_index1,
        int ts_index2,
        int merge_index);
    double compute_weighted_normalized_sum(
        const std::shared_ptr<FactoredTransitionSystem> fts,
        int ts_index1,
        int ts_index2) const;
    void clear();
    void dump_weights() const;
};

class MergeDynamicWeighted : public MergeStrategy {
    const int max_states; // bisimulation option
    const bool use_lr;
    Features *features;
    TaskProxy *task_proxy;

    virtual void dump_strategy_specific_options() const override;
public:
    MergeDynamicWeighted(const Options opts);
    virtual ~MergeDynamicWeighted();
    virtual void initialize(const std::shared_ptr<AbstractTask> task) override;
    virtual std::pair<int, int> get_next(
        std::shared_ptr<FactoredTransitionSystem> fts) override;
    virtual std::string name() const override;
};

#endif
