#include "merge_dynamic_weighted.h"

#include "distances.h"
#include "factored_transition_system.h"
#include "label_equivalence_relation.h"
#include "labels.h"
#include "shrink_bisimulation.h"
#include "transition_system.h"

#include "../causal_graph.h"
#include "../globals.h"
#include "../task_proxy.h"

#include "../options/options.h"

#include "../utils/math.h"

#include <limits>

using namespace std;
using utils::ExitCode;

namespace merge_and_shrink {
const int MINUSINF = numeric_limits<int>::min();

// Helper methods to deal with transition systems

// TODO: copied from DFP
bool is_goal_relevant(const TransitionSystem &ts) {
    int num_states = ts.get_size();
    for (int state = 0; state < num_states; ++state) {
        if (!ts.is_goal_state(state)) {
            return true;
        }
    }
    return false;
}

int compute_number_of_product_transitions(
    const TransitionSystem &ts1, const TransitionSystem &ts2) {
    // TODO: this is copied from the merge constructor of TransitionSystem
    /*
      Note that this computes the number of tranistions in the product
      without considering possible shrinking due to unreachable or
      irrelevant states, which hence may reduce the actual number of
      transitions in the product.
    */
    int number_of_transitions = 0;
    for (const GroupAndTransitions &gat : ts1) {
        const LabelGroup &group1 = gat.label_group;
        const vector<Transition> &transitions1 = gat.transitions;

        // Distribute the labels of this group among the "buckets"
        // corresponding to the groups of ts2.
        unordered_map<int, vector<int>> buckets;
        for (int label_no : group1) {
            int group2_id = ts2.get_group_id_for_label(label_no);
            buckets[group2_id].push_back(label_no);
        }

        // Now create the new groups together with their transitions.
        for (const auto &bucket : buckets) {
            const vector<Transition> &transitions2 =
                ts2.get_transitions_for_group_id(bucket.first);
            int new_transitions_for_new_group = transitions1.size() * transitions2.size();
            number_of_transitions += new_transitions_for_new_group;
        }
    }
    return number_of_transitions;
}

double compute_average_h_value(const Distances &distances) {
    int num_states = distances.get_num_states();
    int sum_distances = 0;
    for (int state = 0; state < num_states; ++state) {
        sum_distances += distances.get_goal_distance(state);
    }
    if (num_states == 0) {
        // For unsolvable transition systems
        return INF;
    }
    return static_cast<double>(sum_distances) / static_cast<double>(num_states);
}

// TODO: exact copy from MergeDFP
void compute_label_ranks(const FactoredTransitionSystem &fts,
                         int index,
                         vector<int> &label_ranks) {
    const TransitionSystem &ts = fts.get_ts(index);
    const Distances &distances = fts.get_dist(index);
    int num_labels = fts.get_num_labels();
    // Irrelevant (and inactive, i.e. reduced) labels have a dummy rank of -1
    label_ranks.resize(num_labels, -1);

    for (const GroupAndTransitions &gat : ts) {
        const LabelGroup &label_group = gat.label_group;
        const vector<Transition> &transitions = gat.transitions;
        // Relevant labels with no transitions have a rank of infinity.
        int label_rank = INF;
        bool group_relevant = false;
        if (static_cast<int>(transitions.size()) == ts.get_size()) {
            /*
              A label group is irrelevant in the earlier notion if it has
              exactly a self loop transition for every state.
            */
            for (const Transition &transition : transitions) {
                if (transition.target != transition.src) {
                    group_relevant = true;
                    break;
                }
            }
        } else {
            group_relevant = true;
        }
        if (!group_relevant) {
            label_rank = -1;
        } else {
            for (const Transition &transition : transitions) {
                label_rank = min(label_rank,
                                 distances.get_goal_distance(transition.target));
            }
        }
        for (int label_no : label_group) {
            label_ranks[label_no] = label_rank;
        }
    }
}

void compute_irrelevant_labels(const FactoredTransitionSystem &fts,
                               vector<vector<bool>> &ts_index_to_irrelevant_labels) {
    int num_ts = fts.get_size();
    ts_index_to_irrelevant_labels.resize(num_ts, vector<bool>());
    int num_labels = fts.get_labels().get_size();
    for (int ts_index = 0; ts_index < num_ts; ++ts_index) {
        if (fts.is_active(ts_index)) {
            vector<bool> irrelevant_labels(num_labels, false);
            const TransitionSystem &ts = fts.get_ts(ts_index);
            for (const GroupAndTransitions &gat : ts) {
                const vector<Transition> &transitions = gat.transitions;
                bool group_relevant = false;
                if (static_cast<int>(transitions.size()) == ts.get_size()) {
                    /*
                      A label group is irrelevant in the earlier notion if it has
                      exactly a self loop transition for every state.
                    */
                    for (size_t i = 0; i < transitions.size(); ++i) {
                        if (transitions[i].target != transitions[i].src) {
                            group_relevant = true;
                            break;
                        }
                    }
                } else {
                    group_relevant = true;
                }
                if (!group_relevant) {
                    const LabelGroup &label_group = gat.label_group;
                    for (int label_no : label_group) {
                        irrelevant_labels[label_no] = true;
                    }
                }
            }
            ts_index_to_irrelevant_labels[ts_index] = irrelevant_labels;
        }
    }
}

// ========================= FEATURE CLASSES ===============================

Feature::Feature(int id, int weight, string name,
                 bool merge_required, bool minimize_value)
    : id(id),
      weight(weight),
      name(name),
      merge_required(merge_required),
      minimize_value(minimize_value) {
}

double Feature::compute_unnormalized_value(
    const FactoredTransitionSystem &fts,
    int ts_index1,
    int ts_index2,
    int merge_index) {
    if (weight) {
        // TODO: get rid of this? we currently also perform this check
        // outside of this method before calling it
        return compute_value(fts, ts_index1, ts_index2, merge_index);
    }
    return 0;
}

void Feature::dump() const {
    cout << name << ": " << weight << endl;
}

CausalConnectionFeature::CausalConnectionFeature(int id, int weight)
    : Feature(id, weight, "causally connected variables", false, false) {
}

double CausalConnectionFeature::compute_value(
    const FactoredTransitionSystem &fts,
    int ts_index1,
    int ts_index2,
    int) {
    // return value in [0,infinity]
    const vector<int> ts1_var_nos = fts.get_ts(ts_index1).get_incorporated_variables();
    const vector<int> ts2_var_nos = fts.get_ts(ts_index2).get_incorporated_variables();
    int edge_count = 0;
    for (int ts1_var_no : ts1_var_nos) {
        for (int ts2_var_no : ts2_var_nos) {
            edge_count += var_pair_causal_links[ts1_var_no][ts2_var_no];
        }
    }
    double max_possible_edges = ts1_var_nos.size() * ts2_var_nos.size() * 2;
    assert(max_possible_edges);
    return static_cast<double>(edge_count) / max_possible_edges;
}

void CausalConnectionFeature::initialize(const TaskProxy &task_proxy, bool dump) {
    const CausalGraph &cg = task_proxy.get_causal_graph();
    int num_variables = task_proxy.get_variables().size();
    var_pair_causal_links.resize(num_variables, vector<int>(num_variables, 0));
    for (VariableProxy var : task_proxy.get_variables()) {
        int var_no = var.get_id();
        const vector<int> &successors = cg.get_successors(var_no);
        for (int succ : successors) {
            var_pair_causal_links[var_no][succ] += 1;
            var_pair_causal_links[succ][var_no] += 1;
        }
        const vector<int> &predecessors = cg.get_predecessors(var_no);
        for (int pred : predecessors) {
            var_pair_causal_links[var_no][pred] += 1;
            var_pair_causal_links[pred][var_no] += 1;
        }
    }
    if (dump) {
        for (int var_no1 = 0; var_no1 < num_variables; ++var_no1) {
            for (int var_no2 = var_no1 + 1; var_no2 < num_variables; ++var_no2) {
                cout << var_no1 << " and " << var_no2 << " have "
                     << var_pair_causal_links[var_no1][var_no2]
                     << "causal graph edges" << endl;
            }
        }
    }
}

BoolCausalConnectionFeature::BoolCausalConnectionFeature(int id, int weight)
    : Feature(id, weight, "boolean causally connected variables", false, false) {
}

double BoolCausalConnectionFeature::compute_value(
    const FactoredTransitionSystem &fts,
    int ts_index1,
    int ts_index2,
    int) {
    // return value in [0,infinity]
    const vector<int> ts1_var_nos = fts.get_ts(ts_index1).get_incorporated_variables();
    const vector<int> ts2_var_nos = fts.get_ts(ts_index2).get_incorporated_variables();
    bool result = false;
    for (int ts1_var_no : ts1_var_nos) {
        for (int ts2_var_no : ts2_var_nos) {
            if (causally_connected_var_pairs[ts1_var_no][ts2_var_no]) {
                result = true;
                break;
            }
        }
    }
    return result;
}

void BoolCausalConnectionFeature::initialize(const TaskProxy &task_proxy, bool dump) {
    const CausalGraph &cg = task_proxy.get_causal_graph();
    int num_variables = task_proxy.get_variables().size();
    causally_connected_var_pairs.resize(num_variables, vector<bool>(num_variables, false));
    for (VariableProxy var : task_proxy.get_variables()) {
        int var_no = var.get_id();
        const vector<int> &successors = cg.get_successors(var_no);
        for (int succ : successors) {
            causally_connected_var_pairs[var_no][succ] = true;
            causally_connected_var_pairs[succ][var_no] = true;
        }
        const vector<int> &predecessors = cg.get_predecessors(var_no);
        for (int pred : predecessors) {
            causally_connected_var_pairs[var_no][pred] = true;
            causally_connected_var_pairs[pred][var_no] = true;
        }
    }
    if (dump) {
        for (int var_no1 = 0; var_no1 < num_variables; ++var_no1) {
            for (int var_no2 = var_no1 + 1; var_no2 < num_variables; ++var_no2) {
                cout << var_no1 << " and " << var_no2 << ": "
                     << (causally_connected_var_pairs[var_no1][var_no2] ? "" : "not ")
                     << "causally connected" << endl;
            }
        }
    }
}

NonAdditivityFeature::NonAdditivityFeature(int id, int weight)
    : Feature(id, weight, "non additive variables", false, false) {
}

double NonAdditivityFeature::compute_value(
    const FactoredTransitionSystem &fts,
    int ts_index1,
    int ts_index2,
    int) {
    // return value in [0,infinity)
    const vector<int> ts1_var_nos = fts.get_ts(ts_index1).get_incorporated_variables();
    const vector<int> ts2_var_nos = fts.get_ts(ts_index2).get_incorporated_variables();
    int not_additive_pair_count = 0;
    for (int ts1_var_no : ts1_var_nos) {
        for (int ts2_var_no : ts2_var_nos) {
            if (!additive_var_pairs[ts1_var_no][ts2_var_no]) {
                ++not_additive_pair_count;
            }
        }
    }
    // NOTE: in contrast to the causally connected variables feature,
    // we consider every pair only once.
    double total_pair_count = ts1_var_nos.size() * ts2_var_nos.size();
    assert(total_pair_count);
    return static_cast<double>(not_additive_pair_count) / total_pair_count;
}

void NonAdditivityFeature::initialize(const TaskProxy &task_proxy, bool dump) {
    int num_variables = task_proxy.get_variables().size();
    additive_var_pairs.resize(num_variables, vector<bool>(num_variables, true));
    for (OperatorProxy op : task_proxy.get_operators()) {
        for (EffectProxy e1 : op.get_effects()) {
            for (EffectProxy e2 : op.get_effects()) {
                int e1_var_id = e1.get_fact().get_variable().get_id();
                int e2_var_id = e2.get_fact().get_variable().get_id();
                additive_var_pairs[e1_var_id][e2_var_id] = false;
                additive_var_pairs[e2_var_id][e1_var_id] = false;
            }
        }
    }
    if (dump) {
        for (int var_no1 = 0; var_no1 < num_variables; ++var_no1) {
            for (int var_no2 = var_no1 + 1; var_no2 < num_variables; ++var_no2) {
                cout << var_no1 << " and " << var_no2 << ": "
                     << (additive_var_pairs[var_no1][var_no2] ? "" : "not ")
                     << "additive" << endl;
            }
        }
    }
}

SmallTransStatesQuotFeature::SmallTransStatesQuotFeature(int id, int weight)
    : Feature(id, weight, "small transitions per states quotient", false, true) {
}

double SmallTransStatesQuotFeature::compute_value(
    const FactoredTransitionSystem &fts,
    int ts_index1,
    int ts_index2,
    int) {
    // return value in [0,infinity)
    const TransitionSystem &ts1 = fts.get_ts(ts_index1);
    const TransitionSystem &ts2 = fts.get_ts(ts_index2);
    double product_states = ts1.get_size() * ts2.get_size();
    double product_transitions = compute_number_of_product_transitions(ts1, ts2);
    if (!product_states) {
        return INF;
    }
    return product_transitions / product_states;
}

HighTransStatesQuotFeature::HighTransStatesQuotFeature(int id, int weight)
    : Feature(id, weight, "high transitions per states quotient", false, false) {
}

double HighTransStatesQuotFeature::compute_value(
    const FactoredTransitionSystem &fts,
    int ts_index1,
    int ts_index2,
    int) {
    // return value in [0,infinity)
    const TransitionSystem &ts1 = fts.get_ts(ts_index1);
    const TransitionSystem &ts2 = fts.get_ts(ts_index2);
    double product_states = ts1.get_size() * ts2.get_size();
    double product_transitions = compute_number_of_product_transitions(ts1, ts2);
    if (!product_states) {
        return INF;
    }
    return product_transitions / product_states;
}

InitHImprovementFeature::InitHImprovementFeature(int id, int weight)
    : Feature(id, weight, "initial h value improvement", true, false) {
}

double InitHImprovementFeature::compute_value(
    const FactoredTransitionSystem &fts,
    int ts_index1,
    int ts_index2,
    int merge_index) {
    // return value in [0,infinity)
    assert(merge_index != -1);
    int new_init_h;
    if (fts.get_ts(merge_index).is_solvable()) {
        new_init_h = fts.get_init_state_goal_distance(merge_index);
    } else {
        // initial state has been pruned
        new_init_h = INF;
    }
    int old_init_h = max(fts.get_init_state_goal_distance(ts_index1),
                         fts.get_init_state_goal_distance(ts_index2));
    int difference = new_init_h - old_init_h;
    if (!difference) {
        return 0;
    }
    if (!old_init_h) {
        return INF;
    }
    return static_cast<double>(difference) / static_cast<double>(old_init_h);
}

AbsoluteInitHFeature::AbsoluteInitHFeature(int id, int weight)
    : Feature(id, weight, "initial h value absolute", true, false) {
}

double AbsoluteInitHFeature::compute_value(
    const FactoredTransitionSystem &fts,
    int,
    int,
    int merge_index) {
    // return value in [0,infinity)
    assert(merge_index != -1);
    if (fts.get_ts(merge_index).is_solvable()) {
        return fts.get_init_state_goal_distance(merge_index);
    } else {
        // initial state has been pruned
        return INF;
    }
}

AbsoluteMaxFFeature::AbsoluteMaxFFeature(int id, int weight)
    : Feature(id, weight, "absolute max f value", true, false) {
}

double AbsoluteMaxFFeature::compute_value(
    const FactoredTransitionSystem &fts,
    int,
    int,
    int merge_index) {
    // return value in [0,infinity)
    assert(merge_index != -1);
    if (fts.get_ts(merge_index).is_solvable()) {
        return fts.get_dist(merge_index).get_max_f();
    } else {
        // initial state has been pruned
        return INF;
    }
}

AbsoluteMaxGFeature::AbsoluteMaxGFeature(int id, int weight)
    : Feature(id, weight, "absolute max g value", true, false) {
}

double AbsoluteMaxGFeature::compute_value(
    const FactoredTransitionSystem &fts,
    int,
    int,
    int merge_index) {
    // return value in [0,infinity)
    assert(merge_index != -1);
    if (fts.get_ts(merge_index).is_solvable()) {
        return fts.get_dist(merge_index).get_max_g();
    } else {
        // initial state has been pruned
        return INF;
    }
}

AbsoluteMaxHFeature::AbsoluteMaxHFeature(int id, int weight)
    : Feature(id, weight, "absolute max h value", true, false) {
}

double AbsoluteMaxHFeature::compute_value(
    const FactoredTransitionSystem &fts,
    int,
    int,
    int merge_index) {
    // return value in [0,infinity)
    assert(merge_index != -1);
    if (fts.get_ts(merge_index).is_solvable()) {
        return fts.get_dist(merge_index).get_max_h();
    } else {
        // initial state has been pruned
        return INF;
    }
}

AvgHImprovementFeature::AvgHImprovementFeature(int id, int weight)
    : Feature(id, weight, "average h value improvement", true, false) {
}

double AvgHImprovementFeature::compute_value(
    const FactoredTransitionSystem &fts,
    int ts_index1,
    int ts_index2,
    int merge_index) {
    // return value in [0,infinity)
    assert(merge_index != -1);
    double new_average_h = compute_average_h_value(fts.get_dist(merge_index));
    double old_average_h = max(compute_average_h_value(fts.get_dist(ts_index1)),
                               compute_average_h_value(fts.get_dist(ts_index2)));
    double difference = new_average_h - old_average_h;
    if (!difference) {
        return 0;
    }
    if (!old_average_h) {
        return INF;
    }
    return static_cast<double>(difference) / static_cast<double>(old_average_h);
}

InitHSumFeature::InitHSumFeature(int id, int weight)
    : Feature(id, weight, "initial h value sum", false, false) {
}

double InitHSumFeature::compute_value(
    const FactoredTransitionSystem &fts,
    int ts_index1,
    int ts_index2,
    int) {
    // return value in [0,infinity)
    int init_h_sum = fts.get_init_state_goal_distance(ts_index1) +
                     fts.get_init_state_goal_distance(ts_index2);
    return init_h_sum;
}

AvgHSumFeature::AvgHSumFeature(int id, int weight)
    : Feature(id, weight, "average h value sum", false, false) {
}

double AvgHSumFeature::compute_value(
    const FactoredTransitionSystem &fts,
    int ts_index1,
    int ts_index2,
    int) {
    // return value in [0,infinity)
    double average_h_sum = compute_average_h_value(fts.get_dist(ts_index1)) +
                           compute_average_h_value(fts.get_dist(ts_index2));
    return average_h_sum;
}

DFPFeature::DFPFeature(int id, int weight)
    : Feature(id, weight, "dfp", false, true) {
}

void DFPFeature::initialize(const TaskProxy &task_proxy, bool) {
    ts_index_to_label_ranks.reserve(task_proxy.get_variables().size() * 2 - 1);
}

//void DFPFeature::precompute_data(const std::FactoredTransitionSystem &fts) override {
//}

void DFPFeature::clear() {
    ts_index_to_label_ranks.clear();
}

double DFPFeature::compute_value(
    const FactoredTransitionSystem &fts,
    int ts_index1,
    int ts_index2,
    int) {
    // return value in [0,infinity)
    if (ts_index_to_label_ranks.empty()) {
        ts_index_to_label_ranks.assign(fts.get_size(), vector<int>());
    }
    const TransitionSystem &ts1 = fts.get_ts(ts_index1);
    const TransitionSystem &ts2 = fts.get_ts(ts_index2);
    int pair_weight = INF;
    if (is_goal_relevant(ts1) || is_goal_relevant(ts2)) {
        vector<int> &label_ranks1 = ts_index_to_label_ranks[ts_index1];
        if (label_ranks1.empty()) {
            compute_label_ranks(fts, ts_index1, label_ranks1);
        }
        vector<int> &label_ranks2 = ts_index_to_label_ranks[ts_index2];
        if (label_ranks2.empty()) {
            compute_label_ranks(fts, ts_index2, label_ranks2);
        }
        assert(!label_ranks1.empty());
        assert(!label_ranks2.empty());
        assert(label_ranks1.size() == label_ranks2.size());
        for (size_t i = 0; i < label_ranks1.size(); ++i) {
            if (label_ranks1[i] != -1 && label_ranks2[i] != -1) {
                // label is relevant in both transition_systems
                int max_label_rank = max(label_ranks1[i], label_ranks2[i]);
                pair_weight = min(pair_weight, max_label_rank);
            }
        }
    }
    // NOTE: if all pairs have infinite weight, the real DFP strategy would
    // prefer merging the "first" pair with at least one goal relevant
    // transition system. This feature here cannot do so.
    return pair_weight;
}

GoalRelevanceFeature::GoalRelevanceFeature(int id, int weight)
    : Feature(id, weight, "goal relevance", false, false) {
}

double GoalRelevanceFeature::compute_value(
    const FactoredTransitionSystem &fts,
    int ts_index1,
    int ts_index2,
    int) {
    // return value in [0,2]
    int pair_weight = 0;
    if (is_goal_relevant(fts.get_ts(ts_index1))) {
        ++pair_weight;
    }
    if (is_goal_relevant(fts.get_ts(ts_index2))) {
        ++pair_weight;
    }
    return pair_weight;
}

NumVariablesFeature::NumVariablesFeature(int id, int weight)
    : Feature(id, weight, "high number of incorporated variables", false, false) {
}

double NumVariablesFeature::compute_value(
    const FactoredTransitionSystem &fts,
    int ts_index1,
    int ts_index2,
    int) {
    // return value in [2,num_variables-1]
    return fts.get_ts(ts_index1).get_incorporated_variables().size() +
           fts.get_ts(ts_index2).get_incorporated_variables().size();
}

ShrinkPerfectlyFeature::ShrinkPerfectlyFeature(int id, int weight)
    : Feature(id, weight, "shrink perfectly", true, false) {
}

double ShrinkPerfectlyFeature::compute_value(
    const FactoredTransitionSystem &fts,
    int,
    int,
    int merge_index) {
    // return value in [0,infinity)
    assert(merge_index != -1);
    if (fts.get_ts(merge_index).is_solvable()) {
        options::Options options;
        options.set<bool>("greedy", false);
        options.set<int>("at_limit", 0);
        ShrinkBisimulation shrink_bisim(options);
        int size_before = fts.get_ts(merge_index).get_size();
        int size_after = shrink_bisim.compute_size_after_perfect_shrink(fts, merge_index);
        assert(size_after <= size_before);
        int difference = size_before - size_after;
        return static_cast<double>(difference) /
               static_cast<double>(size_before);
    } else {
        return INF;
    }
}

NumTransitionsFeature::NumTransitionsFeature(int id, int weight)
    : Feature(id, weight, "small number of transitions", false, true) {
}

double NumTransitionsFeature::compute_value(
    const FactoredTransitionSystem &fts,
    int ts_index1,
    int ts_index2,
    int) {
    // return value in [0,infinity[
    return compute_number_of_product_transitions(fts.get_ts(ts_index1), fts.get_ts(ts_index2));
}

LROpportunitiesFeatures::LROpportunitiesFeatures(int id, int weight)
    : Feature(id, weight, "high number of reducible labels", false, false) {
}

void LROpportunitiesFeatures::clear() {
    ts_pair_to_combinable_label_count.clear();
}

void LROpportunitiesFeatures::precompute_data(
    const FactoredTransitionSystem &fts) {
    // Precompute the set of irrelevant labels for every transition system
    vector<vector<bool>> ts_index_to_irrelevant_labels;
    compute_irrelevant_labels(fts, ts_index_to_irrelevant_labels);

    // Compute the number of labels that are irrelevant in all other transition
    // systems than then current considered pair.
    int num_ts = fts.get_size();
    int num_labels = fts.get_labels().get_size();
    for (int ts_index1 = 0; ts_index1 < num_ts; ++ts_index1) {
        if (fts.is_active(ts_index1)) {
            for (int ts_index2 = ts_index1 + 1; ts_index2 < num_ts; ++ts_index2) {
                if (fts.is_active(ts_index2)) {
                    int count_combinable_labels = 0;
                    for (int label_no = 0; label_no < num_labels; ++label_no) {
                        bool label_irrelevant_in_all_other_ts = true;
                        for (int ts_index3 = 0; ts_index3 < num_ts; ++ts_index3) {
                            if (ts_index3 == ts_index1 || ts_index3 == ts_index2) {
                                continue;
                            }
                            if (fts.is_active(ts_index3)) {
                                if (!ts_index_to_irrelevant_labels[ts_index3][label_no]) {
                                    label_irrelevant_in_all_other_ts = false;
                                    break;
                                }
                            }
                        }
                        if (label_irrelevant_in_all_other_ts) {
                            ++count_combinable_labels;
                        }
                    }
                    ts_pair_to_combinable_label_count[make_pair(ts_index1, ts_index2)] =
                        count_combinable_labels;
                }
            }
        }
    }
}

double LROpportunitiesFeatures::compute_value(
    const FactoredTransitionSystem &,
    int ts_index1,
    int ts_index2,
    int) {
    // return value in [0,infinity[
    int combinable_labels = ts_pair_to_combinable_label_count[make_pair(ts_index1, ts_index2)];
    if (combinable_labels >= 2) {
        // need at least two labels to profit from label reduction
        return combinable_labels;
    }
    return 0;
}

MoreLROpportunitiesFeatures::MoreLROpportunitiesFeatures(int id, int weight)
    : Feature(id, weight, "even higher number of reducible labels", false, false) {
}

void MoreLROpportunitiesFeatures::clear() {
    ts_pair_to_combinable_label_count.clear();
}

void MoreLROpportunitiesFeatures::precompute_data(
    const FactoredTransitionSystem &fts) {
    // Precompute the set of irrelevant labels for every transition system
    vector<vector<bool>> ts_index_to_irrelevant_labels;
    compute_irrelevant_labels(fts, ts_index_to_irrelevant_labels);

    // Compute the number of labels that are locally equivalent in all other
    // transition systems than then current considered pair.
    int num_ts = fts.get_size();
    const Labels &labels = fts.get_labels();
    int num_labels = labels.get_size();
    for (int ts_index1 = 0; ts_index1 < num_ts; ++ts_index1) {
        if (fts.is_active(ts_index1)) {
            for (int ts_index2 = ts_index1 + 1; ts_index2 < num_ts; ++ts_index2) {
                if (fts.is_active(ts_index2)) {
                    int count_combinable_label_pairs = 0;
                    for (int label_no1 = 0; label_no1 < num_labels; ++label_no1) {
                        if (labels.is_current_label(label_no1)) {
                            for (int label_no2 = label_no1 + 1;
                                 label_no2 < num_labels; ++label_no2) {
                                if (labels.is_current_label(label_no2)) {
                                    bool equivalent_in_all_other_ts = true;
                                    for (int ts_index3 = 0; ts_index3 < num_ts; ++ts_index3) {
                                        if (ts_index3 == ts_index1 || ts_index3 == ts_index2) {
                                            continue;
                                        }
                                        if (fts.is_active(ts_index3)) {
                                            const TransitionSystem &ts3 = fts.get_ts(ts_index3);
                                            if (ts3.get_group_id_for_label(label_no1)
                                                != ts3.get_group_id_for_label(label_no2)) {
                                                equivalent_in_all_other_ts = false;
                                                break;
                                            }
                                        }
                                    }
                                    if (equivalent_in_all_other_ts) {
                                        ++count_combinable_label_pairs;
                                    }
                                }
                            }
                        }
                    }
                    ts_pair_to_combinable_label_count[make_pair(ts_index1, ts_index2)] =
                        count_combinable_label_pairs;
                }
            }
        }
    }
}

double MoreLROpportunitiesFeatures::compute_value(
    const FactoredTransitionSystem &,
    int ts_index1,
    int ts_index2,
    int) {
    // return value in [0,infinity[
    int combinable_label_pairs = ts_pair_to_combinable_label_count[make_pair(ts_index1, ts_index2)];
    if (combinable_label_pairs >= 1) {
        // need at least one label pair to profit from label reduction
        return combinable_label_pairs;
    }
    return 0;
}

MIASMFeature::MIASMFeature(int id, int weight)
    : Feature(id, weight, "high number of unreachable and irrelevant states", true, true) {
}

double MIASMFeature::compute_value(
    const FactoredTransitionSystem &fts,
    int ts_index1,
    int ts_index2,
    int merge_index) {
    // return value in [0,infinity)
    assert(merge_index != -1);
    if (fts.get_ts(merge_index).is_solvable()) {
        int expected_size = fts.get_ts(ts_index1).get_size() * fts.get_ts(ts_index2).get_size();
        assert(expected_size);
        int new_size = fts.get_ts(merge_index).get_size();
        assert(new_size <= expected_size);
        return static_cast<double>(new_size) / static_cast<double>(expected_size);
    } else {
        // initial state has been pruned
        // return 0 because this feature is minimized
        return 0;
    }
}

MutexFeature::MutexFeature(int id, int weight)
    : Feature(id, weight, "prefer merging variables that have mutex values ", false, true) {
}

double MutexFeature::compute_value(
    const FactoredTransitionSystem &fts,
    int ts_index1,
    int ts_index2,
    int) {
    // return value in [0,infinity)
    const vector<int> ts1_var_nos = fts.get_ts(ts_index1).get_incorporated_variables();
    const vector<int> ts2_var_nos = fts.get_ts(ts_index2).get_incorporated_variables();
    int mutex_pair_count = 0;
    for (int ts1_var_no : ts1_var_nos) {
        for (int ts2_var_no : ts2_var_nos) {
            if (g_mutex_var_pairs[ts1_var_no][ts2_var_no]) {
                ++mutex_pair_count;
            }
        }
    }
    double total_pair_count = ts1_var_nos.size() * ts2_var_nos.size();
    assert(total_pair_count);
    return static_cast<double>(mutex_pair_count) / total_pair_count;
}

// ========================= FEATURES ====================================

Features::Features(const options::Options opts)
    : debug(opts.get<bool>("debug")),
      merge_required(false) {
    int id = 0;
    features.push_back(new CausalConnectionFeature(
                           id++, opts.get<int>("w_causally_connected_vars")));
    features.push_back(new BoolCausalConnectionFeature(
                           id++, opts.get<int>("w_bool_causally_connected_vars")));
    features.push_back(new NonAdditivityFeature(
                           id++, opts.get<int>("w_nonadditive_vars")));
    features.push_back(new SmallTransStatesQuotFeature(
                           id++, opts.get<int>("w_small_transitions_states_quotient")));
    features.push_back(new HighTransStatesQuotFeature(
                           id++, opts.get<int>("w_high_transitions_states_quotient")));
    features.push_back(new InitHImprovementFeature(
                           id++, opts.get<int>("w_high_initial_h_value_improvement")));
    features.push_back(new AbsoluteInitHFeature(
                           id++, opts.get<int>("w_high_absolute_initial_h_value")));
    features.push_back(new AbsoluteMaxFFeature(
                           id++, opts.get<int>("w_high_absolute_max_f_value")));
    features.push_back(new AbsoluteMaxGFeature(
                           id++, opts.get<int>("w_high_absolute_max_g_value")));
    features.push_back(new AbsoluteMaxHFeature(
                           id++, opts.get<int>("w_high_absolute_max_h_value")));
    features.push_back(new AvgHImprovementFeature(
                           id++, opts.get<int>("w_high_average_h_value_improvement")));
    features.push_back(new InitHSumFeature(
                           id++, opts.get<int>("w_high_initial_h_value_sum")));
    features.push_back(new AvgHSumFeature(
                           id++, opts.get<int>("w_high_average_h_value_sum")));
    features.push_back(new DFPFeature(
                           id++, opts.get<int>("w_dfp")));
    features.push_back(new GoalRelevanceFeature(
                           id++, opts.get<int>("w_goal_relevance")));
    features.push_back(new NumVariablesFeature(
                           id++, opts.get<int>("w_num_variables")));
    features.push_back(new ShrinkPerfectlyFeature(
                           id++, opts.get<int>("w_shrink_perfectly")));
    features.push_back(new NumTransitionsFeature(
                           id++, opts.get<int>("w_num_trans")));
    features.push_back(new LROpportunitiesFeatures(
                           id++, opts.get<int>("w_lr_opp")));
    features.push_back(new MoreLROpportunitiesFeatures(
                           id++, opts.get<int>("w_more_lr_opp")));
    features.push_back(new MIASMFeature(
                           id++, opts.get<int>("w_miasm")));
    features.push_back(new MutexFeature(
                           id++, opts.get<int>("w_mutex")));
    for (Feature *feature : features) {
        if (feature->get_weight() && feature->requires_merge()) {
            merge_required = true;
            break;
        }
    }
}

Features::~Features() {
    for (Feature *feature : features) {
        delete feature;
    }
}

void Features::initialize(const TaskProxy &task_proxy) {
    for (Feature *feature : features) {
        if (feature->get_weight()) {
            feature->initialize(task_proxy, debug);
        }
    }
    clear();
}

void Features::precompute_data(
    const FactoredTransitionSystem &fts) {
    for (Feature *feature : features) {
        if (feature->get_weight()) {
            feature->precompute_data(fts);
        }
    }
}

void Features::update_min_max(int feature_id, double value) {
    // Update minimum an maximum finite values
    if (value == INF || value == MINUSINF) {
        return;
    }
    if (value > max_values[feature_id]) {
        max_values[feature_id] = value;
    }
    if (value < min_values[feature_id]) {
        min_values[feature_id] = value;
    }
}

double Features::normalize_value(int feature_id, double value) const {
    double min = min_values[feature_id];
    double max = max_values[feature_id];
    // Again, only consider finite values. Clamp -infinity and infinity to 0 and 1
    if (value == MINUSINF) {
        return 0;
    }
    if (value == INF) {
        return 1;
    }
    if (max - min == 0) {
        // all three values are the same
        assert(min == value);
        assert(max == value);
        return 0;
    }
    double result = (value - min) / (max - min);
    assert(result >= 0);
    assert(result <= 1);
    return result;
}

void Features::precompute_unnormalized_values(
    const FactoredTransitionSystem &fts,
    int ts_index1,
    int ts_index2,
    int merge_index) {
    vector<double> values;
    values.reserve(features.size());
    for (Feature *feature : features) {
        if (feature->get_weight()) {
            double value = feature->compute_unnormalized_value(fts, ts_index1,
                                                               ts_index2, merge_index);
            update_min_max(feature->get_id(), value);
            values.push_back(value);
        } else {
            // dummy value for correct indices
            // (will be multiplied with weight 0 anyway)
            values.push_back(0);
        }
    }
    unnormalized_values[make_pair(ts_index1, ts_index2)] = values;
}

double Features::compute_weighted_normalized_sum(
    const FactoredTransitionSystem &fts,
    int ts_index1,
    int ts_index2) const {
    pair<int, int> next_pair = make_pair(ts_index1, ts_index2);
    if (!unnormalized_values.count(next_pair)) {
        next_pair = make_pair(ts_index2, ts_index1);
    }
    assert(unnormalized_values.count(next_pair));
    const vector<double> &values = unnormalized_values.at(next_pair);
    double weighted_sum = 0;
    if (debug) {
        cout << "computing weighted normalized sum for "
             << fts.get_ts(ts_index1).tag() << fts.get_ts(ts_index2).tag() << endl;
    }
    for (Feature *feature : features) {
        if (feature->get_weight()) {
            double normalized_value = normalize_value(feature->get_id(),
                                                      values[feature->get_id()]);
            if (feature->minimize()) {
                normalized_value = 1.0 - normalized_value;
            }
            if (debug) {
                cout << "normalized value for feature " << feature->get_name()
                     << ": " << normalized_value << endl;
            }
            weighted_sum += feature->get_weight() * normalized_value;
        }
    }
    if (debug) {
        cout << "weighted normalized sum: " << weighted_sum << endl;
    }
    return weighted_sum;
}

void Features::clear() {
    min_values.assign(features.size(), INF);
    max_values.assign(features.size(), MINUSINF);
    unnormalized_values.clear();
    for (Feature *feature : features) {
        feature->clear();
    }
}

void Features::dump_weights() const {
    for (Feature *feature : features) {
        feature->dump();
    }
}

// ========================= MERGE STRATEGY ===============================

MergeDynamicWeighted::MergeDynamicWeighted(
    FactoredTransitionSystem &fts,
    std::unique_ptr<Features> features,
    std::vector<int> transition_system_order,
    const int max_states,
    const bool use_lr)
    : MergeStrategy(fts),
      features(move(features)),
      transition_system_order(move(transition_system_order)),
      max_states(max_states),
      use_lr(use_lr) {
    if (use_lr) {
        cerr << "Currently not implemented!" << endl;
        utils::exit_with(ExitCode::CRITICAL_ERROR);
    }
}

// TODO: copied from MergeAndShrinkHeuristic
pair<int, int> MergeDynamicWeighted::compute_shrink_sizes(
    int size1, int size2) const {
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

pair<int, int> MergeDynamicWeighted::get_next() {
    int next_index1 = -1;
    int next_index2 = -1;

    features->precompute_data(fts);
    // Go through all transitition systems and compute unnormalized feature values.
    int num_ts = fts.get_size();
    for (int ts_index1 = 0; ts_index1 < num_ts; ++ts_index1) {
        if (fts.is_active(ts_index1)) {
            for (int ts_index2 = ts_index1 + 1; ts_index2 < num_ts; ++ts_index2) {
                if (fts.is_active(ts_index2)) {
                    int copy_ts_index1;
                    int copy_ts_index2;
                    int merge_index = -1;
                    if (features->require_merge()) {
                        // Output for parser
//                            cout << "trying to compute the merge..." << endl;
                        copy_ts_index1 = fts.copy(ts_index1);
                        copy_ts_index2 = fts.copy(ts_index2);
                        pair<int, int> shrink_sizes =
                            compute_shrink_sizes(fts.get_ts(copy_ts_index1).get_size(),
                                                 fts.get_ts(copy_ts_index2).get_size());

                        // shrink before merge (with implicit threshold = 1,
                        // i.e. always try to shrink)
                        options::Options options;
                        options.set<bool>("greedy", false);
                        options.set<int>("at_limit", 0);
                        bool silent = true;
                        ShrinkBisimulation shrink_bisim(options);
                        shrink_bisim.shrink(fts, copy_ts_index1, shrink_sizes.first, silent);
                        shrink_bisim.shrink(fts, copy_ts_index2, shrink_sizes.second, silent);

                        merge_index = fts.merge(copy_ts_index1, copy_ts_index2, true, false);
                        // Output for parser
//                            cout << "...done computing the merge." << endl;
                    }
                    features->precompute_unnormalized_values(fts, ts_index1,
                                                             ts_index2, merge_index);
                    if (features->require_merge()) {
                        // delete and reset
                        fts.release_copies();
                    }
                }
            }
        }
    }

    // Precompute the sorted set of active transition systems
    // TODO: code duplication from MergeDFP again
    assert(!transition_system_order.empty());
    vector<int> sorted_active_ts_indices;
    for (size_t tso_index = 0; tso_index < transition_system_order.size(); ++tso_index) {
        int ts_index = transition_system_order[tso_index];
        if (fts.is_active(ts_index)) {
            sorted_active_ts_indices.push_back(ts_index);
        }
    }

    // Go through all transition systems again and normalize feature values.
    int max_weight = -1;
    unordered_map<int, int> weight_to_count;
    for (size_t i = 0; i < sorted_active_ts_indices.size(); ++i) {
        int ts_index1 = sorted_active_ts_indices[i];
        assert(fts.is_active(ts_index1));
        for (size_t j = i + 1; j < sorted_active_ts_indices.size(); ++j) {
            int ts_index2 = sorted_active_ts_indices[j];
            assert(fts.is_active(ts_index2));
            int pair_weight =
                features->compute_weighted_normalized_sum(fts,
                                                          ts_index1,
                                                          ts_index2);
            if (!weight_to_count.count(pair_weight)) {
                weight_to_count[pair_weight] = 1;
            } else {
                weight_to_count[pair_weight] += 1;
            }
            if (pair_weight > max_weight) {
                max_weight = pair_weight;
                next_index1 = ts_index1;
                next_index2 = ts_index2;
            }
        }
    }
    assert(max_weight != -1);
    /*
      TODO: cache results for all transition systems? only two of the
      transition systems disappear each merge, and only one new arises.
      However, we would need to remove all pairs that the merge transition
      systems were part of, and compute all the pairs with the new one.

      Also, we would need to make sure that at the time that the
      merge strategy is asked for the next pair, the cached results are
      correct, i.e. the transition systems cannot have been shrunk in
      the meantime.
    */
    features->clear();

    int maximum_weight_pair_count = weight_to_count[max_weight];
    assert(maximum_weight_pair_count >= 1);
    if (maximum_weight_pair_count > 1) {
        ++iterations_with_tiebreaking;
        total_tiebreaking_pair_count += maximum_weight_pair_count;
    }

    assert(next_index1 != -1);
    assert(next_index2 != -1);

    return make_pair(next_index1, next_index2);
}
}
