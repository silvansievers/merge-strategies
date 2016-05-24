#include "merge_predefined.h"

#include "factored_transition_system.h"

#include "../option_parser.h"
#include "../plugin.h"

#include "../utils/logging.h"
#include "../utils/system.h"

#include <cassert>
#include <iostream>

using namespace std;
using utils::ExitCode;

namespace merge_and_shrink {
MergePredefined::MergePredefined(
    FactoredTransitionSystem &fts,
    std::vector<std::vector<int>> merge_order)
    : MergeStrategy(fts),
      merge_order(move(merge_order)) {
}

pair<int, int> MergePredefined::get_next() {
    assert(!merge_order.empty());
    int next_index1 = -1;
    int next_index2 = -1;
    const vector<int> &next_pair = merge_order.front();
    assert(next_pair.size() == 2);
    next_index1 = next_pair[0];
    next_index2 = next_pair[1];
    merge_order.erase(merge_order.begin());
    assert(fts.is_active(next_index1));
    assert(fts.is_active(next_index2));
    return make_pair(next_index1, next_index2);
}
}
