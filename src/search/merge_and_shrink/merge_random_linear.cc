#include "merge_random_linear.h"

#include "factored_transition_system.h"

#include "../utils/collections.h"
#include "../utils/rng.h"

#include <cassert>
#include <iostream>

using namespace std;

namespace merge_and_shrink {
MergeRandomLinear::MergeRandomLinear(
    FactoredTransitionSystem &fts,
    vector<int> &&randomized_variable_order)
    : MergeStrategy(fts),
      randomized_variable_order(move(randomized_variable_order)),
      need_first_index(true) {
}

pair<int, int> MergeRandomLinear::get_next() {
    int next_index1;
    if (need_first_index) {
        need_first_index = false;
        next_index1 = randomized_variable_order.front();
        randomized_variable_order.erase(randomized_variable_order.begin());
        cout << "First variable: " << next_index1 << endl;
    } else {
        // The most recent composite transition system is appended at the end
        int num_transition_systems = fts.get_size();
        next_index1 = num_transition_systems - 1;
    }
    int next_index2 = randomized_variable_order.front();
    randomized_variable_order.erase(randomized_variable_order.begin());
    cout << "Next variable: " << next_index2 << endl;
    assert(fts.is_active(next_index1));
    assert(fts.is_active(next_index2));
    return make_pair(next_index1, next_index2);
}
}
