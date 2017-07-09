#ifndef MERGE_AND_SHRINK_UTILS_H
#define MERGE_AND_SHRINK_UTILS_H

#include "types.h"

#include <memory>
#include <vector>

namespace merge_and_shrink {
class Distances;
class FactoredTransitionSystem;
class ShrinkStrategy;
class TransitionSystem;

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
  This method checks if the transition system of the factor at index violates
  the size limit given via new_size (e.g. as computed by compute_shrink_sizes)
  or the threshold shrink_threshold_before_merge that triggers shrinking even
  if the size limit is not violated. If so, trigger the shrinking process.
  Return true iff the factor was actually shrunk.
*/
extern bool shrink_factor(
    FactoredTransitionSystem &fts,
    int index,
    int new_size,
    int shrink_threshold_before_merge,
    const ShrinkStrategy &shrink_strategy,
    Verbosity verbosity);

extern StateEquivalenceRelation compute_pruning_equivalence_relation(
    const TransitionSystem &ts,
    const Distances &distances,
    bool prune_unreachable_states,
    bool prune_irrelevant_states,
    Verbosity verbosity);

/*
  Prune unreachable and/or irrelevant states of the factor at index. This
  requires that init and/or goal distances have been computed accordingly.
  Return true iff any states have been pruned.

  TODO: maybe this functionality belongs to a new class PruningStrategy.
*/
extern bool prune_factor(
    FactoredTransitionSystem &fts,
    int index,
    bool prune_unreachable_states,
    bool prune_irrelevant_states,
    Verbosity verbosity);

/*
  Compute the abstraction mapping based on the given state equivalence
  relation.
*/
extern std::vector<int> compute_abstraction_mapping(
    int num_states,
    const StateEquivalenceRelation &equivalence_relation);

extern bool is_goal_relevant(const TransitionSystem &ts);

extern void shrink_factor(
    const ShrinkStrategy &shrink_strategy,
    TransitionSystem &ts,
    const Distances &dist,
    int new_size,
    Verbosity verbosity);

extern std::pair<std::unique_ptr<TransitionSystem>, std::unique_ptr<Distances>> shrink_and_merge_temporarily(
    const FactoredTransitionSystem &fts,
    int index1,
    int index2,
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
