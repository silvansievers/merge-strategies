#include "utils.h"

#include "factored_transition_system.h"
#include "shrink_bisimulation.h"
#include "transition_system.h"

#include "../options/options.h"

#include "../utils/math.h"

#include <cassert>
#include <math.h>

using namespace std;

namespace merge_and_shrink {
bool is_goal_relevant(const TransitionSystem &ts) {
    int num_states = ts.get_size();
    for (int state = 0; state < num_states; ++state) {
        if (!ts.is_goal_state(state)) {
            return true;
        }
    }
    return false;
}

// TODO: copied from MergeAndShrinkHeuristic
pair<int, int> compute_shrink_sizes(int size1, int size2, int max_states) {
    // Bound both sizes by max allowed size before merge.
    int new_size1 = min(size1, max_states);
    int new_size2 = min(size2, max_states);

    if (!utils::is_product_within_limit(new_size1, new_size2, max_states)) {
        int balanced_size = int(sqrt(max_states));

        if (new_size1 <= balanced_size) {
            // Size of the first transition system is small enough. Use whatever
            // is left for the second transition system.
            new_size2 = max_states / new_size1;
        } else if (new_size2 <= balanced_size) {
            // Inverted case as before.
            new_size1 = max_states / new_size2;
        } else {
            // Both transition systems are too big. We set both target sizes
            // to balanced_size. An alternative would be to set one to
            // N1 = balanced_size and the other to N2 = max_states /
            // balanced_size, to get closer to the allowed maximum.
            // However, this would make little difference (N2 would
            // always be N1, N1 + 1 or N1 + 2), and our solution has the
            // advantage of treating the transition systems symmetrically.
            new_size1 = balanced_size;
            new_size2 = balanced_size;
        }
    }
    assert(new_size1 <= size1 && new_size2 <= size2);
    assert(new_size1 <= max_states);
    assert(new_size2 <= max_states);
    assert(new_size1 * new_size2 <= max_states);
    return make_pair(new_size1, new_size2);
}

int shrink_and_merge_temporarily(
    FactoredTransitionSystem &fts,
    int ts_index1,
    int ts_index2,
    int max_states) {
    // Copy the transition systems (distances etc)
    int copy_ts_index1 = fts.copy(ts_index1);
    int copy_ts_index2 = fts.copy(ts_index2);
    pair<int, int> shrink_sizes =
        compute_shrink_sizes(fts.get_ts(copy_ts_index1).get_size(),
                             fts.get_ts(copy_ts_index2).get_size(),
                             max_states);

    // shrink before merge (with implicit threshold = 1,
    // i.e. always try to shrink)
    options::Options options;
    options.set<bool>("greedy", false);
    options.set<int>("at_limit", 0);
    bool silent = true;
    ShrinkBisimulation shrink_bisim(options);
    shrink_bisim.shrink(fts, copy_ts_index1, shrink_sizes.first, silent);
    shrink_bisim.shrink(fts, copy_ts_index2, shrink_sizes.second, silent);

    // perform the merge
    int merge_index = fts.merge(copy_ts_index1, copy_ts_index2, true, false);
    return merge_index;
}
}
