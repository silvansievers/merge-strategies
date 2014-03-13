#ifndef MERGE_AND_SHRINK_MERGE_SYMMETRIES_H
#define MERGE_AND_SHRINK_MERGE_SYMMETRIES_H

#include "merge_dfp.h"

class Labels;

class MergeSymmetries : public MergeDFP {
protected:
    virtual void dump_strategy_specific_options() const {}
public:
    explicit MergeSymmetries();
    virtual ~MergeSymmetries() {}

    virtual bool done() const;
    virtual std::pair<int, int> get_next(const Labels *labels,
                                         const std::vector<Abstraction *> &all_abstractions);
    virtual std::string name() const;
};

#endif
