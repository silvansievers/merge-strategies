#include "merge_random.h"

#include "factored_transition_system.h"

#include "../utils/collections.h"
#include "../utils/rng.h"

#include <cassert>
#include <iostream>

using namespace std;

namespace merge_and_shrink {
MergeRandom::MergeRandom(
    FactoredTransitionSystem &fts,
    shared_ptr<utils::RandomNumberGenerator> rng)
    : MergeStrategy(fts), rng(move(rng)) {
}

pair<int, int> MergeRandom::get_next() {
    int number_ts = fts.get_size();
    vector<int> active_count_to_ts_index;
    for (int ts_index = 0; ts_index < number_ts; ++ts_index) {
        if (fts.is_active(ts_index)) {
            active_count_to_ts_index.push_back(ts_index);
        }
    }

    utils::RandomNumberGenerator &rng_ = *rng;
    int number_active_ts = active_count_to_ts_index.size();
    int active_index1 = rng_(number_active_ts);
    int active_index2 = rng_(number_active_ts);
    while (active_index2 == active_index1) {
        active_index2 = rng_(number_active_ts);
    }

    assert(utils::in_bounds(active_index1, active_count_to_ts_index));
    assert(utils::in_bounds(active_index2, active_count_to_ts_index));
    int next_index1 = active_count_to_ts_index[active_index1];
    int next_index2 = active_count_to_ts_index[active_index2];

    return make_pair(next_index1, next_index2);
}
}
