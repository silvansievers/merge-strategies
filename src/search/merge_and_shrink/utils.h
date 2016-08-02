#ifndef MERGE_AND_SHRINK_UTILS_H
#define MERGE_AND_SHRINK_UTILS_H

#include <utility>

namespace merge_and_shrink {
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
}

#endif
