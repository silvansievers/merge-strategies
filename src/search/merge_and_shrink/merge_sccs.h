#ifndef MERGE_AND_SHRINK_MERGE_SCCS_H
#define MERGE_AND_SHRINK_MERGE_SCCS_H

#include "merge_dfp.h"

#include "../variable_order_finder.h"

#include <unordered_set>

class Options;

class MergeSCCs : public MergeDFP {
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

    // cg_sccs contain the sccs in order to be merged, from last to first.
    std::vector<std::unordered_set<int> > cg_sccs;
    int number_of_merges_for_scc;
    std::vector<TransitionSystem *> current_transition_systems;
    std::vector<std::pair<int, int> > linear_order;

    std::pair<int, int> get_next_dfp();
protected:
    virtual void dump_strategy_specific_options() const override {}
public:
    MergeSCCs(const Options &options);
    virtual ~MergeSCCs() override = default;
    virtual void initialize(const std::shared_ptr<AbstractTask> task) override;

    virtual std::pair<int, int> get_next(const std::vector<TransitionSystem *> &all_transition_systems) override;
    virtual std::string name() const override;
};

#endif
