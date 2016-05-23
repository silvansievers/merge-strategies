#include "merge_miasm.h"

#include <cassert>

using namespace std;

namespace merge_and_shrink {
MergeMiasm::MergeMiasm(
    FactoredTransitionSystem &fts,
    vector<pair<int, int>> merge_order)
    : MergeStrategy(fts),
      merge_order(move(merge_order)) {
}

pair<int, int> MergeMiasm::get_next() {
    assert(!merge_order.empty());
    pair<int, int> next_merge = merge_order.front();
    merge_order.erase(merge_order.begin());
    return next_merge;
}
}
