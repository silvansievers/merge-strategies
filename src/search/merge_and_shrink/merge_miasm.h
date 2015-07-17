#ifndef MERGE_AND_SHRINK_MERGE_MIASM_H
#define MERGE_AND_SHRINK_MERGE_MIASM_H

#include "merge_miasm_parameters.h"
#include "merge_strategy.h"
//#include "merge_tree.h"
//#include "../operator_cost.h"

#include <cstdlib>

#include <vector>
#include <set>
#include <queue>

class Options;
class TransitionSystem;
class Labels;
class MutexControl;
class SinkSetSearch;

/**
 * @brief The MIASM merging strategy
 * \nosubgrouping
 */
class MergeMiasm : public MergeStrategy {
public:
    /** @brief The option-based constructor */
    MergeMiasm(const Options &opts);
    virtual ~MergeMiasm();
protected:
    /** @name Protected: Options */
    //@{
    /** @brief The option for sink set searching */
    SinkSetSearch *const sink_set_search;
    /** @brief The enum option that specifies the internal merging strategy */
    const MiasmInternal miasm_internal;
    /** @brief The enum option that specifies the external merging strategy */
    const MiasmExternal miasm_external;
    //@}
protected:
    /** @name Protected: Local Computation  */
    //@{
    /** @brief The counter of how many merges have been done */
    std::size_t merge_count;
    //@}
protected:
    /** @name Protected: Result Data */
    //@{
    /** @brief the sink sets */
    std::vector<std::set<int> > sink_sets;
    /** @brief The mutually disjoint subsets which together minimize the
     * ratio of R&R states */
    std::vector<std::set<int> > max_packing;
    /** @brief The list of pairs of indices of intermediate abstractions
     * in the order they are merged in */
    std::vector<std::pair<int, int> > miasm_next;
    //@}
protected:
    virtual void dump_strategy_specific_options() const;
public:
    virtual std::string name() const;
    virtual std::pair<int, int> get_next(
        const std::vector<TransitionSystem *> &all_transition_systems);
    virtual void initialize(const std::shared_ptr<AbstractTask> task) override;
public:
    /**
     * The function that builds a merge-and-shrink abstraction
     * on a subset of variables. It use the merging strategy,
     * the shrinking strategy and label reduction specified
     * for the current MIASM merging strategy
     * @param ordered_varset The subset of variables the abstraction will be
     * build on
     * @param intermediate All abstractions in the final merge tree whose root
     * is the abstraction on the subset
     */
    void build_partial_transition_system(
        const std::vector<int> &ordered_varset,
        std::vector<TransitionSystem *> &intermediate) const;
    /**
     * The greedy method for computing the maximal weighted packing of
     * the family of subsets
     */
    void greedy_max_set_packing();
};

#endif // MERGE_MIASM_H
