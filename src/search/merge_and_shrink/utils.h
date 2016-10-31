#ifndef MERGE_AND_SHRINK_UTILS_H
#define MERGE_AND_SHRINK_UTILS_H

#include <vector>

namespace merge_and_shrink {
class Distances;
class FactoredTransitionSystem;
class ShrinkStrategy;
class TransitionSystem;
enum class Verbosity;

/*
  Compute target sizes for shrinking two transition systems with sizes size1
  and size2 before they are merged. Use the following rules:
  1) Right before merging, the transition systems may have at most
     max_states_before_merge states.
  2) Right after merging, the product has at most max_states_after_merge states.
  3) Transition systems are shrunk as little as necessary to satisfy the above
     constraints. (If possible, neither is shrunk at all.)
  There is often a Pareto frontier of solutions following these rules. In this
  case, balanced solutions (where the target sizes are close to each other)
  are preferred over less balanced ones.
*/
extern std::pair<int, int> compute_shrink_sizes(
    int size1,
    int size2,
    int max_states_before_merge,
    int max_states_after_merge);

/*
  This method checks if the transition system specified via index violates
  the size limit given via new_size (e.g. as computed by compute_shrink_sizes)
  or the threshold shrink_threshold_before_merge that triggers shrinking even
  if the size limit is not violated. If so, the given shrink strategy
  shrink_strategy is used to reduce the size of the transition system to at
  most new_size. Return true iff the transition was modified (i.e. shrunk).
*/
extern bool shrink_transition_system(
    FactoredTransitionSystem &fts,
    int index,
    int new_size,
    int shrink_threshold_before_merge,
    const ShrinkStrategy &shrink_strategy,
    Verbosity verbosity);

extern bool is_goal_relevant(const TransitionSystem &ts);

extern int shrink_and_merge_temporarily(
    FactoredTransitionSystem &fts,
    int ts_index1,
    int ts_index2,
    const ShrinkStrategy &shrink_strategy,
    int max_states,
    int max_states_before_merge,
    int shrink_threshold_before_merge);

extern int compute_number_of_product_transitions(
    const TransitionSystem &ts1, const TransitionSystem &ts2);

extern double compute_average_h_value(const Distances &distances);

extern void compute_irrelevant_labels(
    const FactoredTransitionSystem &fts,
    std::vector<std::vector<bool>> &ts_index_to_irrelevant_labels);
}

#endif
