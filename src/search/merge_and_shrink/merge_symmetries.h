#ifndef MERGE_AND_SHRINK_MERGE_SYMMETRIES_H
#define MERGE_AND_SHRINK_MERGE_SYMMETRIES_H

#include "merge_dfp.h"

#include "../option_parser.h"

class MergeSymmetries : public MergeDFP {
    // options
    const Options options;
    int max_bliss_iterations;

    // statistics
    int iteration_counter;
    int number_of_applied_symmetries;
    bool bliss_limit_reached;
    std::vector<double> bliss_times;

    std::vector<std::pair<int, int> > merge_order; // TODO: change to from last to first?

    void dump_statistics();
protected:
    virtual void dump_strategy_specific_options() const {}
public:
    explicit MergeSymmetries(const Options &options);
    virtual ~MergeSymmetries() {}

    virtual std::pair<int, int> get_next(std::vector<Abstraction *> &all_abstractions);
    virtual std::string name() const;
};

#endif
