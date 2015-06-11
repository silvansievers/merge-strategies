#ifndef MERGE_AND_SHRINK_MERGE_SCCS_H
#define MERGE_AND_SHRINK_MERGE_SCCS_H

#include "merge_dfp.h"

#include <set>

class Options;

class MergeSCCs : public MergeDFP {
    enum SCCOrder {
        TOPOLOGICAL,
        REVERSE_TOPOLOGICAL,
        DECREASING,
        INCREASING
    };
    enum InternalSCCMergeOrder {
        INTERNAL_LINEAR,
        INTERNAL_DFP
    };
    InternalSCCMergeOrder internal_scc_merge_order;
    enum CompositeSCCMergeOrder {
        LINEAR_CREATION_ORDER,
        REVERSE_LINEAR_CREATION_ORDER,
        COMPOSITE_DFP
    };
    CompositeSCCMergeOrder composite_scc_merge_order;

    // TODO: use unordered_set?
    // cg_sccs contain the sccs in order to be merged, from last to first.
    std::vector<std::set<int> > cg_sccs;
    int number_of_merges_for_scc;
    std::vector<TransitionSystem *> current_transition_systems;
    bool merged_sccs;

    std::vector<std::pair<int, int> > internal_linear_order;

    std::pair<int, int> get_next_internal_dfp();
protected:
    virtual void dump_strategy_specific_options() const {}
public:
    MergeSCCs(const Options &options);
    virtual ~MergeSCCs() {}

    virtual std::pair<int, int> get_next(const std::vector<TransitionSystem *> &all_transition_systems);
    virtual std::string name() const;
};

#endif
