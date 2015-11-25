#ifndef MERGE_AND_SHRINK_MERGE_DFP_H
#define MERGE_AND_SHRINK_MERGE_DFP_H

#include "merge_strategy.h"

class Options;
class OptionParser;

class MergeDFP : public MergeStrategy {
    enum AtomicTSOrder {
        REGULAR,
        INVERSE,
        RANDOM1
    };
    AtomicTSOrder atomic_ts_order;
    enum ProductTSOrder {
        OLD_TO_NEW,
        NEW_TO_OLD,
        RANDOM2
    };
    ProductTSOrder product_ts_order;
    bool atomic_before_product;
    bool randomized_order;

    // Store the "DFP" ordering in which transition systems should be considered.
    std::vector<int> transition_system_order;
    void compute_label_ranks(std::shared_ptr<FactoredTransitionSystem> fts,
                             int index,
                             std::vector<int> &label_ranks) const;
    std::pair<int, int> get_next_dfp(
        std::shared_ptr<FactoredTransitionSystem> fts,
        const std::vector<int> &sorted_active_ts_indices) const;
protected:
    virtual void dump_strategy_specific_options() const override {}
public:
    MergeDFP(const Options &options);
    virtual ~MergeDFP() override = default;
    virtual void initialize(const std::shared_ptr<AbstractTask> task) override;

    virtual std::pair<int, int> get_next(std::shared_ptr<FactoredTransitionSystem> fts) override;
    std::pair<int, int> get_next(std::shared_ptr<FactoredTransitionSystem> fts,
                                 const std::vector<int> &ts_indices);
    virtual std::string name() const override;
    static void add_options_to_parser(OptionParser &parser);
};

#endif
