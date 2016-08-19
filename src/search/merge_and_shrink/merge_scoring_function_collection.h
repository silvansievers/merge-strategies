#ifndef MERGE_AND_SHRINK_MERGE_SCORING_FUNCTION_COLLECTION_H
#define MERGE_AND_SHRINK_MERGE_SCORING_FUNCTION_COLLECTION_H

#include "merge_scoring_function.h"

#include "../utils/hash.h"

#include <unordered_map>

namespace options {
class Options;
}

namespace merge_and_shrink {
class MergeScoringFunctionCausalConnection : public MergeScoringFunction {
    std::vector<std::vector<int>> var_pair_causal_links;
protected:
    virtual std::string name() const override;
public:
    MergeScoringFunctionCausalConnection() = default;
    virtual ~MergeScoringFunctionCausalConnection() override = default;
    virtual std::vector<double> compute_scores(
        FactoredTransitionSystem &fts,
        const std::vector<std::pair<int, int>> &merge_candidates) override;
    virtual void initialize(const std::shared_ptr<AbstractTask> &task) override;
};



class MergeScoringFunctionBooleanCausalConnection : public MergeScoringFunction {
    std::vector<std::vector<bool>> causally_connected_var_pairs;
protected:
    virtual std::string name() const override;
public:
    MergeScoringFunctionBooleanCausalConnection() = default;
    virtual ~MergeScoringFunctionBooleanCausalConnection() override = default;
    virtual std::vector<double> compute_scores(
        FactoredTransitionSystem &fts,
        const std::vector<std::pair<int, int>> &merge_candidates) override;
    virtual void initialize(const std::shared_ptr<AbstractTask> &task) override;
};



class MergeScoringFunctionNonAdditivity : public MergeScoringFunction {
    std::vector<std::vector<bool>> additive_var_pairs;
protected:
    virtual std::string name() const override;
public:
    MergeScoringFunctionNonAdditivity() = default;
    virtual ~MergeScoringFunctionNonAdditivity() override = default;
    virtual std::vector<double> compute_scores(
        FactoredTransitionSystem &fts,
        const std::vector<std::pair<int, int>> &merge_candidates) override;
    virtual void initialize(const std::shared_ptr<AbstractTask> &task) override;
};



class MergeScoringFunctionTransitionsStatesQuotient : public MergeScoringFunction {
    bool prefer_high;
protected:
    virtual std::string name() const override;
public:
    explicit MergeScoringFunctionTransitionsStatesQuotient(
        const options::Options &options);
    virtual ~MergeScoringFunctionTransitionsStatesQuotient() override = default;
    virtual std::vector<double> compute_scores(
        FactoredTransitionSystem &fts,
        const std::vector<std::pair<int, int>> &merge_candidates) override;
};



class MergeScoringFunctionInitH : public MergeScoringFunction {
    enum class InitH {
        IMPROVEMENT,
        ABSOLUTE,
        SUM
    };
    InitH inith;
    int max_states;
protected:
    virtual std::string name() const override;
public:
    explicit MergeScoringFunctionInitH(const options::Options &options);
    virtual ~MergeScoringFunctionInitH() override = default;
    virtual std::vector<double> compute_scores(
        FactoredTransitionSystem &fts,
        const std::vector<std::pair<int, int>> &merge_candidates) override;
};



class MergeScoringFunctionMaxFGH : public MergeScoringFunction {
    enum class FGH {
        F,
        G,
        H
    };
    FGH fgh;
    int max_states;
protected:
    virtual std::string name() const override;
public:
    explicit MergeScoringFunctionMaxFGH(const options::Options &options);
    virtual ~MergeScoringFunctionMaxFGH() override = default;
    virtual std::vector<double> compute_scores(
        FactoredTransitionSystem &fts,
        const std::vector<std::pair<int, int>> &merge_candidates) override;
};



class MergeScoringFunctionAvgH : public MergeScoringFunction {
    enum class AvgH {
        IMPROVEMENT,
        ABSOLUTE,
        SUM
    };
    AvgH avgh;
    int max_states;
protected:
    virtual std::string name() const override;
public:
    explicit MergeScoringFunctionAvgH(const options::Options &options);
    virtual ~MergeScoringFunctionAvgH() override = default;
    virtual std::vector<double> compute_scores(
        FactoredTransitionSystem &fts,
        const std::vector<std::pair<int, int>> &merge_candidates) override;
};



class MergeScoringFunctionGoalRelevanceFine : public MergeScoringFunction {
protected:
    virtual std::string name() const override;
public:
    MergeScoringFunctionGoalRelevanceFine() = default;
    virtual ~MergeScoringFunctionGoalRelevanceFine() override = default;
    virtual std::vector<double> compute_scores(
        FactoredTransitionSystem &fts,
        const std::vector<std::pair<int, int>> &merge_candidates) override;
};



class MergeScoringFunctionNumVariables : public MergeScoringFunction {
protected:
    virtual std::string name() const override;
public:
    MergeScoringFunctionNumVariables() = default;
    virtual ~MergeScoringFunctionNumVariables() override = default;
    virtual std::vector<double> compute_scores(
        FactoredTransitionSystem &fts,
        const std::vector<std::pair<int, int>> &merge_candidates) override;
};



class MergeScoringFunctionShrinkPerfectly : public MergeScoringFunction {
    int max_states;
protected:
    virtual std::string name() const override;
public:
    explicit MergeScoringFunctionShrinkPerfectly(
        const options::Options &options);
    virtual ~MergeScoringFunctionShrinkPerfectly() override = default;
    virtual std::vector<double> compute_scores(
        FactoredTransitionSystem &fts,
        const std::vector<std::pair<int, int>> &merge_candidates) override;
};



class MergeScoringFunctionNumTransitions : public MergeScoringFunction {
protected:
    virtual std::string name() const override;
public:
    MergeScoringFunctionNumTransitions() = default;
    virtual ~MergeScoringFunctionNumTransitions() override = default;
    virtual std::vector<double> compute_scores(
        FactoredTransitionSystem &fts,
        const std::vector<std::pair<int, int>> &merge_candidates) override;
};



class MergeScoringFunctionLROpportunities : public MergeScoringFunction {
    std::unordered_map<std::pair<int, int>, int> ts_pair_to_combinable_label_count;
protected:
    virtual std::string name() const override;
public:
    MergeScoringFunctionLROpportunities() = default;
    virtual ~MergeScoringFunctionLROpportunities() override = default;
    virtual std::vector<double> compute_scores(
        FactoredTransitionSystem &fts,
        const std::vector<std::pair<int, int>> &merge_candidates) override;
};



class MergeScoringFunctionMoreLROpportunities : public MergeScoringFunction {
    std::unordered_map<std::pair<int, int>, int> ts_pair_to_combinable_label_count;
protected:
    virtual std::string name() const override;
public:
    MergeScoringFunctionMoreLROpportunities() = default;
    virtual ~MergeScoringFunctionMoreLROpportunities() override = default;
    virtual std::vector<double> compute_scores(
        FactoredTransitionSystem &fts,
        const std::vector<std::pair<int, int>> &merge_candidates) override;
};



class MergeScoringFunctionMutexes : public MergeScoringFunction {
protected:
    virtual std::string name() const override;
public:
    MergeScoringFunctionMutexes() = default;
    virtual ~MergeScoringFunctionMutexes() override = default;
    virtual std::vector<double> compute_scores(
        FactoredTransitionSystem &fts,
        const std::vector<std::pair<int, int>> &merge_candidates) override;
};
}

#endif
