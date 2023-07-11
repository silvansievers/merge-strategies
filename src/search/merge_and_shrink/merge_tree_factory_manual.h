#ifndef MERGE_AND_SHRINK_MERGE_TREE_FACTORY_MANUAL_H
#define MERGE_AND_SHRINK_MERGE_TREE_FACTORY_MANUAL_H

#include "merge_tree_factory.h"

#include <vector>

namespace plugins {
class Options;
}

namespace merge_and_shrink {
class MergeTreeFactoryManual : public MergeTreeFactory {
    std::vector<std::vector<int>> merge_order_list;
    std::string merge_order_tree_string;
protected:
    virtual std::string name() const override;
    virtual void dump_tree_specific_options(utils::LogProxy &log) const override;
public:
    explicit MergeTreeFactoryManual(const plugins::Options &options);
    virtual ~MergeTreeFactoryManual() override = default;
    virtual std::unique_ptr<MergeTree> compute_merge_tree(
        const TaskProxy &task_proxy) override;

    virtual bool requires_init_distances() const override {
        return false;
    }

    virtual bool requires_goal_distances() const override {
        return false;
    }
};
}

#endif
