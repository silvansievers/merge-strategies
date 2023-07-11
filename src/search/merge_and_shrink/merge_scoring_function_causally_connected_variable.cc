#include "merge_scoring_function_causally_connected_variable.h"

#include "factored_transition_system.h"
#include "transition_system.h"

#include "../task_proxy.h"

#include "../plugins/plugin.h"

#include "../tasks/root_task.h"

#include "../task_utils/causal_graph.h"

using namespace std;

namespace merge_and_shrink {
MergeScoringFunctionCausallyConnectedVariable::MergeScoringFunctionCausallyConnectedVariable()
    : // HACK!
      task_proxy(*tasks::g_root_task) {
}

vector<double> MergeScoringFunctionCausallyConnectedVariable::compute_scores(
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

    vector<int> composite_vars;
    for (int ts_index : composite_indices) {
        const vector<int> &vars = fts.get_transition_system(ts_index).get_incorporated_variables();
        composite_vars.insert(composite_vars.end(), vars.begin(), vars.end());
    }

    int num_vars = task_proxy.get_variables().size();
    vector<bool> is_causal_predecessor(num_vars, false);
    const causal_graph::CausalGraph &cg = task_proxy.get_causal_graph();
    for (int var : composite_vars) {
        const vector<int> &connected_vars = cg.get_eff_to_pre(var);
        for (int connected_var : connected_vars)
            is_causal_predecessor[connected_var] = true;
    }

    vector<double> scores;
    scores.reserve(merge_candidates.size());
    for (pair<int, int> merge_candidate : merge_candidates) {
        int ts_index1 = merge_candidate.first;
        int ts_index2 = merge_candidate.second;
        int score = INF;
        const vector<int> &variables1 =
            fts.get_transition_system(ts_index1).get_incorporated_variables();
        const vector<int> &variables2 =
            fts.get_transition_system(ts_index2).get_incorporated_variables();
        if ((variables1.size() == 1 && is_causal_predecessor[variables1.front()])
            || (variables2.size() == 1 && is_causal_predecessor[variables2.front()])) {
            score = 0;
        }
        scores.push_back(score);
    }
    return scores;
}

void MergeScoringFunctionCausallyConnectedVariable::initialize(
    const TaskProxy &task) {
    initialized = true;
    this->task_proxy = task;
}

string MergeScoringFunctionCausallyConnectedVariable::name() const {
    return "causally connected variable";
}

class MergeScoringFunctionCausallyConnectedVariableFeature : public plugins::TypedFeature<MergeScoringFunction, MergeScoringFunctionCausallyConnectedVariable> {
public:
    MergeScoringFunctionCausallyConnectedVariableFeature() : TypedFeature("causally_connected_variable") {
        document_title("CausallyConnectedVariable");
        document_synopsis(
            "This scoring function assumes that all non-linear merge candidates "
            "have been filtered out before. It assigns a merge candidate a value "
            "of 0 iff at least one of the components is atomic and its variable "
            "causally connected to (at least one of) the composite transition "
            "system(s), and otherwise infinity."
            );
    }

    virtual shared_ptr<MergeScoringFunctionCausallyConnectedVariable> create_component(
        const plugins::Options &, const utils::Context &) const override {
        return make_shared<MergeScoringFunctionCausallyConnectedVariable>();
    }
};

static plugins::FeaturePlugin<MergeScoringFunctionCausallyConnectedVariableFeature> _plugin;
}
