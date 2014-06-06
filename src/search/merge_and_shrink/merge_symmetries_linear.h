#ifndef MERGE_AND_SHRINK_MERGE_SYMMETRIES_LINEAR_H
#define MERGE_AND_SHRINK_MERGE_SYMMETRIES_LINEAR_H

#include "merge_linear.h"

#include "../option_parser.h"

#include <set>

class MergeSymmetriesLinear : public MergeLinear {
    const Options options;
    std::set<int> abs_to_merge;
    bool started_merging_for_symmetries;
    // the following serves for statistics output
    // number of applications of at least one symmetry
    int number_of_applied_symmetries;
    int atomic_symmetries; // symmetries affecting one abstraction
    int binary_symmetries; // symmetries affecting two abstractions
    int other_symmetries; // symmetries affecting more than two abstractions
    int iteration_counter;
    int max_symmetry_iterations;
    void dump_statistics() const;
protected:
    virtual void dump_strategy_specific_options() const {}
public:
    explicit MergeSymmetriesLinear(const Options &options);
    virtual ~MergeSymmetriesLinear() {}

    virtual bool done() const;
    virtual std::pair<int, int> get_next(const std::vector<Abstraction *> &all_abstractions);
    virtual std::string name() const;
};

#endif
