#ifndef MERGE_AND_SHRINK_MERGE_RANDOM_LINEAR_H
#define MERGE_AND_SHRINK_MERGE_RANDOM_LINEAR_H

#include "merge_strategy.h"

namespace options {
class Options;
}
namespace utils {
    class RandomNumberGenerator;
}

namespace merge_and_shrink {
class MergeRandomLinear : public MergeStrategy {
    int random_seed;
    std::shared_ptr<utils::RandomNumberGenerator> rng;
    bool need_first_index;
    std::vector<int> randomized_variable_order;
protected:
    virtual void dump_strategy_specific_options() const override;
public:
    explicit MergeRandomLinear(const options::Options &options);
    virtual ~MergeRandomLinear() override = default;
    virtual void initialize(const std::shared_ptr<AbstractTask> task) override;

    virtual std::pair<int, int> get_next(
        FactoredTransitionSystem &fts) override;
    virtual std::string name() const override;
};
}

#endif
