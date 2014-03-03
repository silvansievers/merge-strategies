#ifndef MERGE_AND_SHRINK_SYMMETRIES_SYMMETRIES_VARIABLE_ORDER_FINDER_H
#define MERGE_AND_SHRINK_SYMMETRIES_SYMMETRIES_VARIABLE_ORDER_FINDER_H

#include "../variable_order_finder.h"

#include "../merge_and_shrink_heuristic.h" // needed for MergeStrategy type;
// TODO: move that type somewhere else?
#include "symmetries.h"

#include <set>
#include <vector>

// TODO: rename class? Here also symmetries are computed
class SymmetriesVariableOrderFinder : public VariableOrderFinder{
private:
    Symmetries symmetries;
    std::set<int> abs_to_merge;
    bool found_symmetry;
    void set_used_var(int var_no);
public:
    SymmetriesVariableOrderFinder(MergeStrategy merge_strategy,
                                  bool is_first,
                                  bool debug_graph_creator,
                                  int version,
                                  const std::vector<Abstraction *> &atomic_abstractions);
    void compute_symmetries(const std::vector<Abstraction *> &abstractions);
    int get_next_to_merge();
    void apply_symmetries(const std::vector<Abstraction *> &abstractions);

    bool found_symmetries() const {return found_symmetry; }
    bool can_apply_now() const {return found_symmetry && abs_to_merge.empty(); }
};

#endif // MERGE_AND_SHRINK_SYMMETRIES_SYMMETRIES_VARIABLE_ORDER_FINDER_H
