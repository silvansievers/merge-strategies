#ifndef MERGE_TREE_H
#define MERGE_TREE_H

#include "merge_miasm_parameters.h"

#include "subset_info.h"

#include "../../ext/tree.hh"

#include <cstdlib>
#include <memory>
#include <vector>
#include <set>

class AbstractTask;
class CausalGraph;

namespace merge_and_shrink {
class ComparatorSortPacking : public ComparatorVarSet {
public:
    ComparatorSortPacking(const std::shared_ptr<AbstractTask> task,
                          const MiasmExternal &ext_,
                          const VarSetInfoRegistry *p_si_ = 0);
    virtual ~ComparatorSortPacking();
    virtual bool operator()(
        const std::set<int> &set_i, const std::set<int> &set_j) const;
protected:
    const MiasmExternal &ext;
};

class MergeTree {
public:
    void get_order(std::vector<std::pair<int, int> > &merge_next_, int num_vars);
    int get_slot(const tree<std::set<int> >::iterator ti);
    tree<std::set<int>> &get_tree() {
        return merge_tree;
    }
protected:
    tree<std::set<int> > merge_tree;
    int slot_count;
    std::vector<std::pair<int, int> > merge_next;
protected:
    static void merge_subs(const tree<std::set<int> > &left,
                           const tree<std::set<int> > &right,
                           tree<std::set<int> > &merged);
};

class MiasmMergeTree : public MergeTree {
public:
    MiasmMergeTree(const std::vector<std::set<int> > &packing_,
                   const MiasmInternal internal_,
                   const MiasmExternal external_,
                   const VarSetInfoRegistry *p_si,
                   const std::shared_ptr<AbstractTask> task);
private:
    std::vector<std::set<int> > packing;
    const MiasmInternal internal;
    const MiasmExternal external;
    const std::shared_ptr<AbstractTask> task;
    const CausalGraph &causal_graph;
private:
    void get_internal_tree(const std::set<int> &varset,
                           tree<std::set<int> > &internal_tree);
private:
    void update_pred(const std::size_t i);
    std::set<int> pred;
    std::set<int> goal;
};
}

#endif // MERGE_TREE_H
