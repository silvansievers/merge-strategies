#include "merge_scoring_function_linear.h"

#include "factored_transition_system.h"

#include "../options/option_parser.h"
#include "../options/plugin.h"

#include <iostream>

using namespace std;

namespace merge_and_shrink {
vector<int> MergeScoringFunctionLinear::compute_scores(
    FactoredTransitionSystem &fts,
    const vector<pair<int, int>> &merge_candidates) {
    vector<int> scores;
    scores.reserve(merge_candidates.size());
    int composite_ts_index = fts.get_size() - 1;
    for (const pair<int, int> &merge : merge_candidates) {
        if (merge.first == composite_ts_index ||
            merge.second == composite_ts_index) {
            scores.push_back(0);
        } else {
            scores.push_back(INF);
        }
    }
    return scores;
}

string MergeScoringFunctionLinear::name() const {
    return "linear";
}

static shared_ptr<MergeScoringFunction>_parse(options::OptionParser &parser) {
    parser.document_synopsis(
        "Linear",
        "This scoring function assigns all pairs containing the newest "
        "composite transition system a score of 0, and ifinity otherwise. "
        "The only exception is before the first merge, where there is no "
        "composite transition system, and hence all candidates receive a "
        "score of infinity.\n"
        "Note that if there are several composite transition systems, the "
        "strategy will not be linear anymore. This may still be useful if "
        "used on a subset of transition systems to be merge linearly.");

    if (parser.dry_run())
        return nullptr;
    else
        return make_shared<MergeScoringFunctionLinear>();
}

static options::PluginShared<MergeScoringFunction> _plugin("linear", _parse);
}
