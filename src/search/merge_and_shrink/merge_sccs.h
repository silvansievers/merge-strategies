#ifndef MERGE_AND_SHRINK_MERGE_SCCS_H
#define MERGE_AND_SHRINK_MERGE_SCCS_H

#include "merge_strategy.h"

#include "../variable_order_finder.h"

#include <unordered_set>

class MergeDFP;
class Options;

class MergeSCCs : public MergeStrategy {
    enum SCCOrder {
        TOPOLOGICAL,
        REVERSE_TOPOLOGICAL,
        DECREASING,
        INCREASING
    };
    SCCOrder scc_order;
    enum MergeOrder {
        LINEAR,
        DFP
    };
    MergeOrder merge_order;
    VariableOrderType var_order_type;
    MergeDFP *merge_dfp;

    // cg_sccs contain the sccs in order to be merged, from last to first.
    std::vector<std::unordered_set<int>> cg_sccs;
    int number_of_merges_for_scc;
    std::vector<int> current_scc_ts_indices;
    std::vector<std::pair<int, int>> linear_order;

    std::pair<int, int> get_next_dfp(
        std::shared_ptr<FactoredTransitionSystem> fts);
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
