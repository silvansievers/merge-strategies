#ifndef MERGE_AND_SHRINK_MERGE_SCCS_H
#define MERGE_AND_SHRINK_MERGE_SCCS_H

#include "merge_dfp.h"

#include <set>

class Options;

class MergeSCCs : public MergeDFP {
    std::vector<std::set<int> > cg_sccs;
    int number_of_merges_for_scc;
    std::vector<TransitionSystem *> current_transition_systems;

    std::pair<int, int> get_next_current_scc();
protected:
    virtual void dump_strategy_specific_options() const {}
public:
    MergeSCCs();
    virtual ~MergeSCCs() {}

    virtual std::pair<int, int> get_next(const std::vector<TransitionSystem *> &all_transition_systems);
    virtual std::string name() const;
};

#endif
