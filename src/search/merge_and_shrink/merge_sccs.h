#ifndef MERGE_AND_SHRINK_MERGE_SCCS_H
#define MERGE_AND_SHRINK_MERGE_SCCS_H

#include "merge_strategy.h"

#include <unordered_set>

class MergeDFP;
class Options;

class MergeSCCs : public MergeStrategy {
    const Options *options;
    enum OrderOfSCCs {
        TOPOLOGICAL,
        REVERSE_TOPOLOGICAL,
        DECREASING,
        INCREASING
    };
    OrderOfSCCs order_of_sccs;
    enum InternalMergeOrder {
        LINEAR1,
        DFP1,
        LINEAR_MANUAL
    };
    InternalMergeOrder internal_merge_order;
    enum MergedSCCsMergeOrder {
        LINEAR2,
        DFP2
    };
    MergedSCCsMergeOrder merged_sccs_merge_order;
    std::vector<int> linear_variable_order;
    MergeDFP *merge_dfp;

    std::vector<bool> is_causally_linked;
    std::vector<bool> is_goal_variable;
    std::vector<std::vector<int>> cg_successors;
    std::vector<std::vector<int>> cg_predecessors;

    // cg_sccs contain the sccs in order to be merged, from last to first.
    std::vector<std::unordered_set<int>> cg_sccs;
    int number_of_merges_for_scc; // TODO: get rid of, use current_scc_ts_indices.size()
    std::vector<int> current_scc_ts_indices;
    bool merged_all_sccs;
    std::vector<int> indices_of_merged_sccs;
    bool start_merging_sccs;
protected:
    virtual void dump_strategy_specific_options() const override {}
public:
    MergeSCCs(const Options &options);
    virtual ~MergeSCCs();
    virtual void initialize(const std::shared_ptr<AbstractTask> task) override;

    virtual std::pair<int, int> get_next(
        std::shared_ptr<FactoredTransitionSystem> fts) override;
    virtual std::string name() const override;
};

#endif
