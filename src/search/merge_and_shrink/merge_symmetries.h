#ifndef MERGE_AND_SHRINK_MERGE_SYMMETRIES_H
#define MERGE_AND_SHRINK_MERGE_SYMMETRIES_H

#include "merge_dfp.h"

#include "../option_parser.h"

class MergeSymmetries : public MergeDFP {
    const Options options;
    //int max_iterations;
    std::vector<std::pair<int, int> > merge_order; // TODO: change to from last to first?

    int number_of_applied_symmetries;
    //int iteration_counter;
    bool bliss_limit_reached;

    void dump_statistics() const;
protected:
    virtual void dump_strategy_specific_options() const {}
public:
    explicit MergeSymmetries(const Options &options);
    virtual ~MergeSymmetries() {}

    virtual std::pair<int, int> get_next(std::vector<Abstraction *> &all_abstractions);
    virtual std::string name() const;
};

#endif
