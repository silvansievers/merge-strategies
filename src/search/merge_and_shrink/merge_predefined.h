#ifndef MERGE_AND_SHRINK_MERGE_PREDEFINED_H
#define MERGE_AND_SHRINK_MERGE_PREDEFINED_H

#include "merge_strategy.h"

namespace options {
class Options;
}

namespace merge_and_shrink {
class BinaryTree;

class MergePredefined : public MergeStrategy {
    std::vector<std::vector<int>> merge_order;
    BinaryTree *root;
protected:
    virtual void dump_strategy_specific_options() const override;
public:
    explicit MergePredefined(const options::Options &options);
    virtual ~MergePredefined() override;
    virtual void initialize(const std::shared_ptr<AbstractTask> task);

    virtual std::pair<int, int> get_next(
        FactoredTransitionSystem &fts) override;
    virtual std::string name() const override;
};
}

#endif