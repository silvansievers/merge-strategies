#ifndef MERGE_AND_SHRINK_MERGE_STRATEGY_FACTORY_MIASM_H
#define MERGE_AND_SHRINK_MERGE_STRATEGY_FACTORY_MIASM_H

#include "merge_strategy_factory.h"

#include "miasm/merge_tree.h"

#include <vector>
#include <set>

namespace options {
class OptionParser;
class Options;
}

namespace merge_and_shrink {
class MiasmMergeTree;

/**
 * @brief The MIASM merging strategy
 * \nosubgrouping
 */
class MergeStrategyFactoryMiasm : public MergeStrategyFactory {
private:
    /**
     * The greedy method for computing the maximal weighted packing of
     * the family of subsets
     */
    void greedy_max_set_packing();
    /** @name Protected: Options */
    options::Options *options;
    //@{
    /** @brief The enum option that specifies the internal merging strategy */
    const MiasmInternal miasm_internal;
    /** @brief The enum option that specifies the external merging strategy */
    const MiasmExternal miasm_external;
    //@}
    /** @name Protected: Local Computation  */
    //@{
    /** @name Protected: Result Data */
    //@{
    /** @brief the sink sets */
    std::vector<std::set<int> > sink_sets;
    /** @brief The mutually disjoint subsets which together minimize the
     * ratio of R&R states */
    std::vector<std::set<int> > max_packing;
    //@}
protected:
    virtual void dump_strategy_specific_options() const override;
public:
    /** @brief The option-based constructor */
    explicit MergeStrategyFactoryMiasm(const options::Options &opts);
    virtual ~MergeStrategyFactoryMiasm() override;
    virtual std::string name() const;
    MiasmMergeTree *compute_merge_tree(const std::shared_ptr<AbstractTask> task);
    virtual std::unique_ptr<MergeStrategy> compute_merge_strategy(
        const std::shared_ptr<AbstractTask> task,
        FactoredTransitionSystem &fts) override;
    static void add_options_to_parser(options::OptionParser &parser);
};
}

#endif
