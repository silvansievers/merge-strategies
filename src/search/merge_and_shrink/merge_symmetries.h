#ifndef MERGE_AND_SHRINK_MERGE_SYMMETRIES_H
#define MERGE_AND_SHRINK_MERGE_SYMMETRIES_H

#include "merge_dfp.h"

#include "../option_parser.h"

#include <set>

class MergeSymmetries : public MergeDFP {
    const Options options;
    std::set<int> abs_to_merge;
    bool started_merging_for_symmetries;
    // the following serves for statistics output
    int atomic_symmetries; // symmetries affecting one abstraction
    int binary_symmetries; // symmetries affecting two abstractions
    int other_symmetries; // symmetries affecting more than two abstractions
    int iteration_counter;
    int max_symmetry_iterations;
    void dump_statistics() const;
protected:
    virtual void dump_strategy_specific_options() const {}
public:
    explicit MergeSymmetries(const Options &options);
    virtual ~MergeSymmetries() {}

    virtual bool done() const;
    virtual std::pair<int, int> get_next(const std::vector<Abstraction *> &all_abstractions);
    virtual std::string name() const;
};

#endif
