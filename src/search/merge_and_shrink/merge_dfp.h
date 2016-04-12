#ifndef MERGE_AND_SHRINK_MERGE_DFP_H
#define MERGE_AND_SHRINK_MERGE_DFP_H

#include "merge_strategy.h"

namespace options {
class OptionParser;
class Options;
}
class TaskProxy;

namespace merge_and_shrink {
enum class AtomicTSOrder {
    REGULAR,
    INVERSE,
    RANDOM
};

enum class ProductTSOrder {
    OLD_TO_NEW,
    NEW_TO_OLD,
    RANDOM
};

class MergeDFP : public MergeStrategy {
    AtomicTSOrder atomic_ts_order;
    ProductTSOrder product_ts_order;
    bool atomic_before_product;
    bool randomized_order;

    // Store the "DFP" ordering in which transition systems should be considered.
    std::vector<int> transition_system_order;
    void compute_label_ranks(const FactoredTransitionSystem &fts,
                             int index,
                             std::vector<int> &label_ranks) const;
    std::pair<int, int> get_next_dfp(
        const FactoredTransitionSystem &fts,
        const std::vector<int> &sorted_active_ts_indices);
protected:
    virtual void dump_strategy_specific_options() const override {}
public:
    MergeDFP(const options::Options &options);
    virtual ~MergeDFP() override = default;
    virtual void initialize(const std::shared_ptr<AbstractTask> task) override;

    virtual std::pair<int, int> get_next(FactoredTransitionSystem &fts) override;
    std::pair<int, int> get_next(FactoredTransitionSystem &fts,
                                 const std::vector<int> &sorted_indices);
    virtual std::string name() const override;
    static void compute_ts_order(const TaskProxy &task_proxy,
                                 AtomicTSOrder atomic_ts_order,
                                 ProductTSOrder product_ts_order,
                                 bool atomic_before_product,
                                 bool randomized_order,
                                 std::vector<int> &ts_order);
    static void add_options_to_parser(options::OptionParser &parser, bool dfp_defaults = true);
};
}

#endif
