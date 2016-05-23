#include "merge_strategy.h"

using namespace std;

namespace merge_and_shrink {
MergeStrategy::MergeStrategy(
    FactoredTransitionSystem &fts)
    : fts(fts), iterations_with_tiebreaking(0), total_tiebreaking_pair_count(0) {
}
}
