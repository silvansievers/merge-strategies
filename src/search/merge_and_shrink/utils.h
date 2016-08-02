#ifndef MERGE_AND_SHRINK_UTILS_H
#define MERGE_AND_SHRINK_UTILS_H

#include <vector>

namespace merge_and_shrink {
class Distances;
class FactoredTransitionSystem;
class TransitionSystem;

extern bool is_goal_relevant(const TransitionSystem &ts);

extern std::pair<int, int> compute_shrink_sizes(
    int size1, int size2, int max_states);

extern int shrink_and_merge_temporarily(
    FactoredTransitionSystem &fts,
    int ts_index1,
    int ts_index2,
    int max_states);

extern int compute_number_of_product_transitions(
    const TransitionSystem &ts1, const TransitionSystem &ts2);

extern double compute_average_h_value(const Distances &distances);

extern void compute_irrelevant_labels(
    const FactoredTransitionSystem &fts,
    std::vector<std::vector<bool>> &ts_index_to_irrelevant_labels);
}

#endif
