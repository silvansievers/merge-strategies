#include "merge_scoring_function_linear.h"

#include "factored_transition_system.h"
#include "transition_system.h"

#include "../plugins/plugin.h"

#include <iostream>
#include <set>

using namespace std;

namespace merge_and_shrink {
vector<double> MergeScoringFunctionLinear::compute_scores(
    const FactoredTransitionSystem &fts,
    const vector<pair<int, int>> &merge_candidates) {
    set<int> composite_indices;
    for (int ts_index = 0; ts_index < fts.get_size(); ++ts_index) {
        if (fts.is_active(ts_index)) {
            int num_vars = fts.get_transition_system(ts_index).get_incorporated_variables().size();
            if (num_vars > 1) {
                composite_indices.insert(ts_index);
            }
        }
    }

    vector<double> scores;
    scores.reserve(merge_candidates.size());
    for (const pair<int, int> &merge : merge_candidates) {
        if (composite_indices.count(merge.first) ||
            composite_indices.count(merge.second)) {
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

class MergeScoringFunctionLinearFeature : public plugins::TypedFeature<MergeScoringFunction, MergeScoringFunctionLinear> {
public:
    MergeScoringFunctionLinearFeature() : TypedFeature("sf_linear") {
        document_title("Linear");
        document_synopsis(
            "This scoring function assigns all pairs containing the newest "
            "composite transition system a score of 0, and ifinity otherwise. "
            "The only exception is before the first merge, where there is no "
            "composite transition system, and hence all candidates receive a "
            "score of infinity.\n"
            "Note that if there are several composite transition systems, the "
            "strategy will not be linear anymore. This may still be useful if "
            "used on a subset of transition systems to be merge linearly.");
    }

    virtual shared_ptr<MergeScoringFunctionLinear> create_component(
        const plugins::Options &, const utils::Context &) const override {
        return make_shared<MergeScoringFunctionLinear>();
    }
};

static plugins::FeaturePlugin<MergeScoringFunctionLinearFeature> _plugin;
}
