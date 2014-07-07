#ifndef MERGE_AND_SHRINK_MERGE_SYMMETRIES_H
#define MERGE_AND_SHRINK_MERGE_SYMMETRIES_H

#include "merge_dfp.h"

#include "../option_parser.h"

class MergeSymmetries : public MergeDFP {
    const Options options;
    int max_iterations;
    enum InternalMerging {
        LINEAR,
        NON_LINEAR
    };
    InternalMerging internal_merging;

    std::vector<int> abs_to_merge;
    bool started_merging_for_symmetries;

    // the following serves for statistics output
    // number of applications of at least one symmetry
    int number_of_applied_symmetries;
    int atomic_symmetries; // symmetries affecting one abstraction
    int binary_symmetries; // symmetries affecting two abstractions
    int other_symmetries; // symmetries affecting more than two abstractions
    int iteration_counter;

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
