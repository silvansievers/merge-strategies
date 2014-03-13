#ifndef MERGE_AND_SHRINK_MERGE_SYMMETRIES_H
#define MERGE_AND_SHRINK_MERGE_SYMMETRIES_H

#include "merge_dfp.h"

#include <set>

class Labels;
class Options;

class MergeSymmetries : public MergeDFP {
    const Options &options;
    class Symmetries *symmetries;
    std::set<int> abs_to_merge;
    /* this is the index at which the composite abstraction we aim to
       construct by merging all indices in abs_to_merge is stored in
       all_abstractions (this is currently the *first* index we return from
       abs_to_merge. */
    int index_of_composite_abs;
protected:
    virtual void dump_strategy_specific_options() const {}
public:
    explicit MergeSymmetries(const Options &options);
    virtual ~MergeSymmetries() {}

    virtual bool done() const;
    virtual std::pair<int, int> get_next(const std::vector<Abstraction *> &all_abstractions);
    virtual std::string name() const;

    void initialize(const Labels *labels);
};

#endif
