#ifndef MERGE_AND_SHRINK_MERGE_MIASM_H
#define MERGE_AND_SHRINK_MERGE_MIASM_H

#include "merge_miasm_parameters.h"

#include "../merge_strategy.h"

#include "../../option_parser.h"

#include <set>

class Options;

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
    const Options options;
    //@{
    /** @brief The enum option that specifies the internal merging strategy */
    const MiasmInternal miasm_internal;
    /** @brief The enum option that specifies the external merging strategy */
    const MiasmExternal miasm_external;
    //@}
    /** @name Protected: Local Computation  */
    //@{
    /** @brief The counter of how many merges have been done */
    std::size_t merge_count;
    //@}
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
    virtual void dump_strategy_specific_options() const;
public:
    virtual std::string name() const;
    virtual std::pair<int, int> get_next(
        std::shared_ptr<FactoredTransitionSystem> fts) override;
    virtual void initialize(const std::shared_ptr<AbstractTask> task) override;
    /**
     * The greedy method for computing the maximal weighted packing of
     * the family of subsets
     */
    void greedy_max_set_packing();
};

#endif // MERGE_MIASM_H
