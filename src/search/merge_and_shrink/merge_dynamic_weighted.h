#ifndef MERGE_AND_SHRINK_MERGE_DYNAMIC_WEIGHTED_H
#define MERGE_AND_SHRINK_MERGE_DYNAMIC_WEIGHTED_H

#include "merge_strategy.h"

#include "../utils/hash.h"

#include <memory>
#include <unordered_map>

namespace options {
class Options;
}
class TaskProxy;

namespace merge_and_shrink {
const int NUM_FEATURES = 7;

class Feature {
    const int id;
    const int weight;
    const std::string name;
    const bool merge_required;
    const bool minimize_value;
    virtual double compute_value(
        const FactoredTransitionSystem &fts,
        int ts_index1,
        int ts_index2,
        int merge_index) = 0;
public:
    Feature(int id, int weight, std::string name,
            bool requires_merge, bool minimize_value);
    virtual ~Feature() {}
    virtual void initialize(const TaskProxy &, bool) {}
    virtual void precompute_data(
        const FactoredTransitionSystem &) {}
    double compute_unnormalized_value(
        const FactoredTransitionSystem &fts,
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
        const FactoredTransitionSystem &fts,
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
        const FactoredTransitionSystem &fts,
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
        const FactoredTransitionSystem &fts,
        int ts_index1,
        int ts_index2,
        int merge_index) override;
public:
    NonAdditivityFeature(int id, int weight);
    virtual void initialize(const TaskProxy &task_proxy, bool dump) override;
};

class SmallTransStatesQuotFeature : public Feature {
    virtual double compute_value(
        const FactoredTransitionSystem &fts,
        int ts_index1,
        int ts_index2,
        int merge_index) override;
public:
    SmallTransStatesQuotFeature(int id, int weight);
};

class HighTransStatesQuotFeature : public Feature {
    virtual double compute_value(
        const FactoredTransitionSystem &fts,
        int ts_index1,
        int ts_index2,
        int merge_index) override;
public:
    HighTransStatesQuotFeature(int id, int weight);
};

class InitHImprovementFeature : public Feature {
    virtual double compute_value(
        const FactoredTransitionSystem &fts,
        int ts_index1,
        int ts_index2,
        int merge_index) override;
public:
    InitHImprovementFeature(int id, int weight);
};

class AbsoluteInitHFeature : public Feature {
    virtual double compute_value(
        const FactoredTransitionSystem &fts,
        int ts_index1,
        int ts_index2,
        int merge_index) override;
public:
    AbsoluteInitHFeature(int id, int weight);
};

class AbsoluteMaxFFeature : public Feature {
    virtual double compute_value(
        const FactoredTransitionSystem &fts,
        int ts_index1,
        int ts_index2,
        int merge_index) override;
public:
    AbsoluteMaxFFeature(int id, int weight);
};

class AbsoluteMaxGFeature : public Feature {
    virtual double compute_value(
        const FactoredTransitionSystem &fts,
        int ts_index1,
        int ts_index2,
        int merge_index) override;
public:
    AbsoluteMaxGFeature(int id, int weight);
};

class AbsoluteMaxHFeature : public Feature {
    virtual double compute_value(
        const FactoredTransitionSystem &fts,
        int ts_index1,
        int ts_index2,
        int merge_index) override;
public:
    AbsoluteMaxHFeature(int id, int weight);
};

class AvgHImprovementFeature : public Feature {
    virtual double compute_value(
        const FactoredTransitionSystem &fts,
        int ts_index1,
        int ts_index2,
        int merge_index) override;
public:
    AvgHImprovementFeature(int id, int weight);
};

class InitHSumFeature : public Feature {
    virtual double compute_value(
        const FactoredTransitionSystem &fts,
        int ts_index1,
        int ts_index2,
        int merge_index) override;
public:
    InitHSumFeature(int id, int weight);
};

class AvgHSumFeature : public Feature {
    virtual double compute_value(
        const FactoredTransitionSystem &fts,
        int ts_index1,
        int ts_index2,
        int merge_index) override;
public:
    AvgHSumFeature(int id, int weight);
};

class GoalRelevanceFeature : public Feature {
    virtual double compute_value(
        const FactoredTransitionSystem &fts,
        int ts_index1,
        int ts_index2,
        int merge_index) override;
public:
    GoalRelevanceFeature(int id, int weight);
};

class NumVariablesFeature : public Feature {
    virtual double compute_value(
        const FactoredTransitionSystem &fts,
        int ts_index1,
        int ts_index2,
        int merge_index) override;
public:
    NumVariablesFeature(int id, int weight);
};

class ShrinkPerfectlyFeature : public Feature {
    virtual double compute_value(
        const FactoredTransitionSystem &fts,
        int ts_index1,
        int ts_index2,
        int merge_index) override;
public:
    ShrinkPerfectlyFeature(int id, int weight);
};

class NumTransitionsFeature : public Feature {
    virtual double compute_value(
        const FactoredTransitionSystem &fts,
        int ts_index1,
        int ts_index2,
        int merge_index) override;
public:
    NumTransitionsFeature(int id, int weight);
};

class LROpportunitiesFeatures : public Feature {
    std::unordered_map<std::pair<int, int>, int> ts_pair_to_combinable_label_count;
    virtual double compute_value(
        const FactoredTransitionSystem &fts,
        int ts_index1,
        int ts_index2,
        int merge_index) override;
public:
    LROpportunitiesFeatures(int id, int weight);
    virtual void clear() override;
    void precompute_data(const FactoredTransitionSystem &fts) override;
};

class MoreLROpportunitiesFeatures : public Feature {
    std::unordered_map<std::pair<int, int>, int> ts_pair_to_combinable_label_count;
    virtual double compute_value(
        const FactoredTransitionSystem &fts,
        int ts_index1,
        int ts_index2,
        int merge_index) override;
public:
    MoreLROpportunitiesFeatures(int id, int weight);
    virtual void clear() override;
    void precompute_data(const FactoredTransitionSystem &fts) override;
};

class MutexFeature : public Feature {
    virtual double compute_value(
        const FactoredTransitionSystem &fts,
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
    explicit Features(const options::Options opts);
    ~Features();
    void initialize(const TaskProxy &task_proxy);
    void precompute_data(const FactoredTransitionSystem &fts);
    bool require_merge() const {
        return merge_required;
    }
    void precompute_unnormalized_values(
        const FactoredTransitionSystem &fts,
        int ts_index1,
        int ts_index2,
        int merge_index);
    double compute_weighted_normalized_sum(
        const FactoredTransitionSystem &fts,
        int ts_index1,
        int ts_index2) const;
    void clear();
    void dump_weights() const;
};

class MergeDynamicWeighted : public MergeStrategy {
    std::unique_ptr<Features> features;
    std::vector<int> transition_system_order; // TS order as in DFP
    const int max_states; // limit to compute shrink sizes
    const bool use_lr;

    int iterations_with_tiebreaking;
    int total_tiebreaking_pair_count;

    std::pair<int, int> compute_shrink_sizes(int size1, int size2) const;
public:
    MergeDynamicWeighted(
        FactoredTransitionSystem &fts,
        std::unique_ptr<Features> features,
        std::vector<int> transition_system_order,
        const int max_states,
        const bool use_lr);
    virtual ~MergeDynamicWeighted() override = default;
    virtual std::pair<int, int> get_next() override;
    virtual std::pair<int, int> get_dfp_tiebreaking_statistics() const override;
};
}

#endif
