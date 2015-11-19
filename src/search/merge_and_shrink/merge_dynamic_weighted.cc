#include "merge_dynamic_weighted.h"

#include "distances.h"
#include "factored_transition_system.h"
#include "labels.h"
#include "shrink_bisimulation.h"
#include "transition_system.h"

#include "../causal_graph.h"
#include "../globals.h"
#include "../option_parser.h"
#include "../plugin.h"
#include "../task_proxy.h"


using namespace std;

const int MINUSINF = numeric_limits<int>::min();

// Helper methods to deal with transition systems

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
    for (TSConstIterator group1_it = ts1.begin();
         group1_it != ts1.end(); ++group1_it) {
        // Distribute the labels of this group among the "buckets"
        // corresponding to the groups of ts2.
        unordered_map<int, vector<int>> buckets;
        for (LabelConstIter label_it = group1_it.begin();
             label_it != group1_it.end(); ++label_it) {
            int label_no = *label_it;
            int group2_id = ts2.get_group_id_for_label(label_no);
            buckets[group2_id].push_back(label_no);
        }
        // Now buckets contains all equivalence classes that are
        // refinements of group1.

        // Now create the new groups together with their transitions.
        const vector<Transition> &transitions1 = group1_it.get_transitions();
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
void compute_label_ranks(shared_ptr<FactoredTransitionSystem> fts,
                         int index,
                         vector<int> &label_ranks)  {
    const TransitionSystem &ts = fts->get_ts(index);
    const Distances &distances = fts->get_dist(index);
    int num_labels = fts->get_num_labels();
    // Irrelevant (and inactive, i.e. reduced) labels have a dummy rank of -1
    label_ranks.resize(num_labels, -1);

    for (TSConstIterator group_it = ts.begin();
         group_it != ts.end(); ++group_it) {
        // Relevant labels with no transitions have a rank of infinity.
        int label_rank = INF;
        const vector<Transition> &transitions = group_it.get_transitions();
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
            label_rank = -1;
        } else {
            for (size_t i = 0; i < transitions.size(); ++i) {
                const Transition &t = transitions[i];
                label_rank = min(label_rank, distances.get_goal_distance(t.target));
            }
        }
        for (LabelConstIter label_it = group_it.begin();
             label_it != group_it.end(); ++label_it) {
            int label_no = *label_it;
            label_ranks[label_no] = label_rank;
        }
    }
}

void compute_irrelevant_labels(const shared_ptr<FactoredTransitionSystem> fts,
                               vector<vector<bool>> &ts_index_to_irrelevant_labels) {
    int num_ts = fts->get_size();
    ts_index_to_irrelevant_labels.resize(num_ts, vector<bool>());
    int num_labels = fts->get_labels()->get_size();
    for (int ts_index = 0; ts_index < num_ts; ++ts_index) {
        if (fts->is_active(ts_index)) {
            vector<bool> irrelevant_labels(num_labels, false);
            const TransitionSystem &ts = fts->get_ts(ts_index);
            for (TSConstIterator group_it = ts.begin();
                 group_it != ts.end(); ++group_it) {
                const vector<Transition> &transitions = group_it.get_transitions();
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
                    for (LabelConstIter label_it = group_it.begin();
                         label_it != group_it.end(); ++label_it) {
                        int label_no = *label_it;
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
    const shared_ptr<FactoredTransitionSystem> fts,
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
    const shared_ptr<FactoredTransitionSystem> fts,
    int ts_index1,
    int ts_index2,
    int) {
    // return value in [0,infinity]
    const vector<int> ts1_var_nos = fts->get_ts(ts_index1).get_incorporated_variables();
    const vector<int> ts2_var_nos = fts->get_ts(ts_index2).get_incorporated_variables();
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
    const shared_ptr<FactoredTransitionSystem> fts,
    int ts_index1,
    int ts_index2,
    int) {
    // return value in [0,infinity]
    const vector<int> ts1_var_nos = fts->get_ts(ts_index1).get_incorporated_variables();
    const vector<int> ts2_var_nos = fts->get_ts(ts_index2).get_incorporated_variables();
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
    const shared_ptr<FactoredTransitionSystem> fts,
    int ts_index1,
    int ts_index2,
    int) {
    // return value in [0,infinity)
    const vector<int> ts1_var_nos = fts->get_ts(ts_index1).get_incorporated_variables();
    const vector<int> ts2_var_nos = fts->get_ts(ts_index2).get_incorporated_variables();
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
    const shared_ptr<FactoredTransitionSystem> fts,
    int ts_index1,
    int ts_index2,
    int) {
    // return value in [0,infinity)
    const TransitionSystem &ts1 = fts->get_ts(ts_index1);
    const TransitionSystem &ts2 = fts->get_ts(ts_index2);
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
    const shared_ptr<FactoredTransitionSystem> fts,
    int ts_index1,
    int ts_index2,
    int) {
    // return value in [0,infinity)
    const TransitionSystem &ts1 = fts->get_ts(ts_index1);
    const TransitionSystem &ts2 = fts->get_ts(ts_index2);
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
    const shared_ptr<FactoredTransitionSystem> fts,
    int ts_index1,
    int ts_index2,
    int merge_index) {
    // return value in [0,infinity)
    assert(merge_index != -1);
    int new_init_h;
    if (fts->get_ts(merge_index).is_solvable()) {
        new_init_h = fts->get_init_state_goal_distance(merge_index);
    } else {
        // initial state has been pruned
        new_init_h = INF;
    }
    int old_init_h = max(fts->get_init_state_goal_distance(ts_index1),
                         fts->get_init_state_goal_distance(ts_index2));
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
    const shared_ptr<FactoredTransitionSystem> fts,
    int,
    int,
    int merge_index) {
    // return value in [0,infinity)
    assert(merge_index != -1);
    if (fts->get_ts(merge_index).is_solvable()) {
        return fts->get_init_state_goal_distance(merge_index);
    } else {
        // initial state has been pruned
        return INF;
    }
}

AbsoluteMaxFFeature::AbsoluteMaxFFeature(int id, int weight)
    : Feature(id, weight, "absolute max f value", true, false) {
}

double AbsoluteMaxFFeature::compute_value(
    const shared_ptr<FactoredTransitionSystem> fts,
    int,
    int,
    int merge_index) {
    // return value in [0,infinity)
    assert(merge_index != -1);
    if (fts->get_ts(merge_index).is_solvable()) {
        return fts->get_dist(merge_index).get_max_f();
    } else {
        // initial state has been pruned
        return INF;
    }
}

AbsoluteMaxGFeature::AbsoluteMaxGFeature(int id, int weight)
    : Feature(id, weight, "absolute max g value", true, false) {
}

double AbsoluteMaxGFeature::compute_value(
    const shared_ptr<FactoredTransitionSystem> fts,
    int,
    int,
    int merge_index) {
    // return value in [0,infinity)
    assert(merge_index != -1);
    if (fts->get_ts(merge_index).is_solvable()) {
        return fts->get_dist(merge_index).get_max_g();
    } else {
        // initial state has been pruned
        return INF;
    }
}

AbsoluteMaxHFeature::AbsoluteMaxHFeature(int id, int weight)
    : Feature(id, weight, "absolute max h value", true, false) {
}

double AbsoluteMaxHFeature::compute_value(
    const shared_ptr<FactoredTransitionSystem> fts,
    int,
    int,
    int merge_index) {
    // return value in [0,infinity)
    assert(merge_index != -1);
    if (fts->get_ts(merge_index).is_solvable()) {
        return fts->get_dist(merge_index).get_max_h();
    } else {
        // initial state has been pruned
        return INF;
    }
}

AvgHImprovementFeature::AvgHImprovementFeature(int id, int weight)
    : Feature(id, weight, "average h value improvement", true, false) {
}

double AvgHImprovementFeature::compute_value(
    const shared_ptr<FactoredTransitionSystem> fts,
    int ts_index1,
    int ts_index2,
    int merge_index) {
    // return value in [0,infinity)
    assert(merge_index != -1);
    double new_average_h = compute_average_h_value(fts->get_dist(merge_index));
    double old_average_h = max(compute_average_h_value(fts->get_dist(ts_index1)),
                               compute_average_h_value(fts->get_dist(ts_index2)));
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
    const shared_ptr<FactoredTransitionSystem> fts,
    int ts_index1,
    int ts_index2,
    int) {
    // return value in [0,infinity)
    int init_h_sum = fts->get_init_state_goal_distance(ts_index1) +
                     fts->get_init_state_goal_distance(ts_index2);
    return init_h_sum;
}

AvgHSumFeature::AvgHSumFeature(int id, int weight)
    : Feature(id, weight, "average h value sum", false, false) {
}

double AvgHSumFeature::compute_value(
    const shared_ptr<FactoredTransitionSystem> fts,
    int ts_index1,
    int ts_index2,
    int) {
    // return value in [0,infinity)
    double average_h_sum = compute_average_h_value(fts->get_dist(ts_index1)) +
                           compute_average_h_value(fts->get_dist(ts_index2));
    return average_h_sum;
}

DFPFeature::DFPFeature(int id, int weight)
    : Feature(id, weight, "dfp", false, true) {
}

void DFPFeature::initialize(const TaskProxy &task_proxy, bool) {
    ts_index_to_label_ranks.reserve(task_proxy.get_variables().size() * 2 - 1);
}

//void DFPFeature::precompute_data(const std::shared_ptr<FactoredTransitionSystem> fts) override {
//}

void DFPFeature::clear() {
    ts_index_to_label_ranks.clear();
}

double DFPFeature::compute_value(
    const shared_ptr<FactoredTransitionSystem> fts,
    int ts_index1,
    int ts_index2,
    int) {
    // return value in [0,infinity)
    if (ts_index_to_label_ranks.empty()) {
        ts_index_to_label_ranks.assign(fts->get_size(), vector<int>());
    }
    const TransitionSystem &ts1 = fts->get_ts(ts_index1);
    const TransitionSystem &ts2 = fts->get_ts(ts_index2);
    int pair_weight = INF;
    if (ts1.is_goal_relevant() || ts2.is_goal_relevant()) {
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
    const shared_ptr<FactoredTransitionSystem> fts,
    int ts_index1,
    int ts_index2,
    int) {
    // return value in [0,2]
    int pair_weight = 0;
    if (fts->get_ts(ts_index1).is_goal_relevant()) {
        ++pair_weight;
    }
    if (fts->get_ts(ts_index2).is_goal_relevant()) {
        ++pair_weight;
    }
    return pair_weight;
}

NumVariablesFeature::NumVariablesFeature(int id, int weight)
    : Feature(id, weight, "high number of incorporated variables", false, false) {
}

double NumVariablesFeature::compute_value(
    const shared_ptr<FactoredTransitionSystem> fts,
    int ts_index1,
    int ts_index2,
    int) {
    // return value in [2,num_variables-1]
    return fts->get_ts(ts_index1).get_incorporated_variables().size() +
           fts->get_ts(ts_index2).get_incorporated_variables().size();
}

ShrinkPerfectlyFeature::ShrinkPerfectlyFeature(int id, int weight)
    : Feature(id, weight, "shrink perfectly", true, false) {
}

double ShrinkPerfectlyFeature::compute_value(
    const shared_ptr<FactoredTransitionSystem> fts,
    int,
    int ,
    int merge_index) {
    // return value in [0,infinity)
    assert(merge_index != -1);
    if (fts->get_ts(merge_index).is_solvable()) {
        Options options;
        options.set<int>("max_states", INF);
        options.set<int>("max_states_before_merge", INF);
        options.set<int>("threshold", 1);
        options.set<bool>("greedy", false);
        options.set<int>("at_limit", 0);
        ShrinkBisimulation shrink_bisim(options);
        int size_before = fts->get_ts(merge_index).get_size();
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
    const shared_ptr<FactoredTransitionSystem> fts,
    int ts_index1,
    int ts_index2,
    int) {
    // return value in [0,infinity[
    return compute_number_of_product_transitions(fts->get_ts(ts_index1), fts->get_ts(ts_index2));
}

LROpportunitiesFeatures::LROpportunitiesFeatures(int id, int weight)
    : Feature(id, weight, "high number of reducible labels", false, false) {
}

void LROpportunitiesFeatures::clear() {
    ts_pair_to_combinable_label_count.clear();
}

void LROpportunitiesFeatures::precompute_data(
    const shared_ptr<FactoredTransitionSystem> fts) {
    // Precompute the set of irrelevant labels for every transition system
    vector<vector<bool>> ts_index_to_irrelevant_labels;
    compute_irrelevant_labels(fts, ts_index_to_irrelevant_labels);

    // Compute the number of labels that are irrelevant in all other transition
    // systems than then current considered pair.
    int num_ts = fts->get_size();
    int num_labels = fts->get_labels()->get_size();
    for (int ts_index1 = 0; ts_index1 < num_ts; ++ts_index1) {
        if (fts->is_active(ts_index1)) {
            for (int ts_index2 = ts_index1 + 1; ts_index2 < num_ts; ++ts_index2) {
                if (fts->is_active(ts_index2)) {
                    int count_combinable_labels = 0;
                    for (int label_no = 0; label_no < num_labels; ++label_no) {
                        bool label_irrelevant_in_all_other_ts = true;
                        for (int ts_index3 = 0; ts_index3 < num_ts; ++ts_index3) {
                            if (ts_index3 == ts_index1 || ts_index3 == ts_index2) {
                                continue;
                            }
                            if (fts->is_active(ts_index3)) {
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
    const shared_ptr<FactoredTransitionSystem>,
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
    const shared_ptr<FactoredTransitionSystem> fts) {
    // Precompute the set of irrelevant labels for every transition system
    vector<vector<bool>> ts_index_to_irrelevant_labels;
    compute_irrelevant_labels(fts, ts_index_to_irrelevant_labels);

    // Compute the number of labels that are locally equivalent in all other
    // transition systems than then current considered pair.
    int num_ts = fts->get_size();
    const shared_ptr<Labels> labels = fts->get_labels();
    int num_labels = labels->get_size();
    for (int ts_index1 = 0; ts_index1 < num_ts; ++ts_index1) {
        if (fts->is_active(ts_index1)) {
            for (int ts_index2 = ts_index1 + 1; ts_index2 < num_ts; ++ts_index2) {
                if (fts->is_active(ts_index2)) {
                    int count_combinable_label_pairs = 0;
                    for (int label_no1 = 0; label_no1 < num_labels; ++label_no1) {
                        if (labels->is_current_label(label_no1)) {
                            for (int label_no2 = label_no1 + 1;
                                 label_no2 < num_labels; ++label_no2) {
                                if (labels->is_current_label(label_no2)) {
                                    bool equivalent_in_all_other_ts = true;
                                    for (int ts_index3 = 0; ts_index3 < num_ts; ++ts_index3) {
                                        if (ts_index3 == ts_index1 || ts_index3 == ts_index2) {
                                            continue;
                                        }
                                        if (fts->is_active(ts_index3)) {
                                            const TransitionSystem &ts3 = fts->get_ts(ts_index3);
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
    const shared_ptr<FactoredTransitionSystem>,
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
    const shared_ptr<FactoredTransitionSystem> fts,
    int ts_index1,
    int ts_index2,
    int merge_index) {
    // return value in [0,infinity)
    assert(merge_index != -1);
    if (fts->get_ts(merge_index).is_solvable()) {
        int expected_size = fts->get_ts(ts_index1).get_size() * fts->get_ts(ts_index2).get_size();
        assert(expected_size);
        int new_size = fts->get_ts(merge_index).get_size();
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
    const shared_ptr<FactoredTransitionSystem> fts,
    int ts_index1,
    int ts_index2,
    int) {
    // return value in [0,infinity)
    const vector<int> ts1_var_nos = fts->get_ts(ts_index1).get_incorporated_variables();
    const vector<int> ts2_var_nos = fts->get_ts(ts_index2).get_incorporated_variables();
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

Features::Features(const Options opts)
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
    const shared_ptr<FactoredTransitionSystem> fts) {
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
    const shared_ptr<FactoredTransitionSystem> fts,
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
    const shared_ptr<FactoredTransitionSystem> fts,
    int ts_index1,
    int ts_index2) const {
    const vector<double> &values = unnormalized_values.at(make_pair(ts_index1, ts_index2));
    double weighted_sum = 0;
    if (debug) {
        cout << "computing weighted normalized sum for "
             << fts->get_ts(ts_index1).tag() << fts->get_ts(ts_index2).tag() << endl;
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

MergeDynamicWeighted::MergeDynamicWeighted(const Options opts)
    : MergeStrategy(),
      max_states(opts.get<int>("max_states")),
      use_lr(opts.get<bool>("use_lr")) {
    if (use_lr) {
        cerr << "Currently not implemented!" << endl;
        exit_with(EXIT_CRITICAL_ERROR);
    }
    features = new Features(opts);
}

MergeDynamicWeighted::~MergeDynamicWeighted() {
    delete task_proxy;
    delete features;
}

void MergeDynamicWeighted::dump_strategy_specific_options() const {
    features->dump_weights();
}

void MergeDynamicWeighted::initialize(const shared_ptr<AbstractTask> task) {
    MergeStrategy::initialize(task);
    task_proxy = new TaskProxy(*task);
    features->initialize(*task_proxy);
}

pair<int, int> MergeDynamicWeighted::get_next(
    shared_ptr<FactoredTransitionSystem> fts) {
    int next_index1 = -1;
    int next_index2 = -1;

    int num_ts = fts->get_size();
    if (remaining_merges == 1) {
        for (int ts_index = 0; ts_index < num_ts - 1; ++ts_index) {
            if (fts->is_active(ts_index)) {
                next_index1 = ts_index;
                break;
            }
        }
        next_index2 = num_ts - 1; // the previously added transition system
        assert(next_index2 != next_index1);
        assert(fts->is_active(next_index2));
    } else {
        features->precompute_data(fts);
        // Go through all transitition systems and compute unnormalized feature values.
        for (int ts_index1 = 0; ts_index1 < num_ts; ++ts_index1) {
            if (fts->is_active(ts_index1)) {
                for (int ts_index2 = ts_index1 + 1; ts_index2 < num_ts; ++ts_index2) {
                    if (fts->is_active(ts_index2)) {
                        int copy_ts_index1;
                        int copy_ts_index2;
                        int merge_index = -1;
                        if (features->require_merge()) {
                            // Output for parser
                            cout << "trying to compute the merge..." << endl;
                            copy_ts_index1 = fts->copy(ts_index1);
                            copy_ts_index2 = fts->copy(ts_index2);
                            Options options;
                            options.set<int>("max_states", max_states);
                            options.set<int>("max_states_before_merge", max_states);
                            options.set<int>("threshold", 1);
                            options.set<bool>("greedy", false);
                            options.set<int>("at_limit", 0);
                            ShrinkBisimulation shrink_bisim(options);
                            shrink_bisim.shrink(fts, copy_ts_index1, copy_ts_index2, true);
                            merge_index = fts->merge(copy_ts_index1, copy_ts_index2, true, false);
                            // Output for parser
                            cout << "...done computing the merge." << endl;
                        }
                        features->precompute_unnormalized_values(fts, ts_index1,
                                                                 ts_index2, merge_index);
                        if (features->require_merge()) {
                            // delete and reset
                            fts->release_copies();
                        }
                    }
                }
            }
        }

        // Go through all transition systems again and normalize feature values.
        int max_weight = -1;
        for (int ts_index1 = 0; ts_index1 < num_ts; ++ts_index1) {
            if (fts->is_active(ts_index1)) {
                for (int ts_index2 = ts_index1 + 1; ts_index2 < num_ts; ++ts_index2) {
                    if (fts->is_active(ts_index2)) {
                        int pair_weight =
                            features->compute_weighted_normalized_sum(fts,
                                                                      ts_index1,
                                                                      ts_index2);
                        if (pair_weight > max_weight) {
                            max_weight = pair_weight;
                            next_index1 = ts_index1;
                            next_index2 = ts_index2;
                        }
                    }
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
    }

    assert(next_index1 != -1);
    assert(next_index2 != -1);

    --remaining_merges;
    return make_pair(next_index1, next_index2);
}

string MergeDynamicWeighted::name() const {
    return "dynamic merging";
}

static shared_ptr<MergeStrategy>_parse(OptionParser &parser) {
    parser.add_option<bool>(
        "debug", "debug", "false");
    parser.add_option<int>("max_states", "shrink strategy option", "50000");
    parser.add_option<bool>("use_lr", "use label reduction", "false");
    parser.add_option<int>(
        "w_causally_connected_vars",
        "prefer merging variables that are causally connected ",
        "0",
        Bounds("0", "100"));
    parser.add_option<int>(
        "w_bool_causally_connected_vars",
        "prefer merging variables that are causally connected ",
        "0",
        Bounds("0", "100"));
    parser.add_option<int>(
        "w_nonadditive_vars",
        "avoid merging additive variables",
        "0",
        Bounds("0", "100"));
    parser.add_option<int>(
        "w_small_transitions_states_quotient",
        "prefer merging 'sparse' transition systems",
        "0",
        Bounds("0", "100"));
    parser.add_option<int>(
        "w_high_transitions_states_quotient",
        "prefer merging 'dense' transition systems",
        "0",
        Bounds("0", "100"));
    parser.add_option<int>(
        "w_high_initial_h_value_improvement",
        "prefer merging transition systems with high initial h value improvement",
        "0",
        Bounds("0", "100"));
    parser.add_option<int>(
        "w_high_absolute_initial_h_value",
        "prefer merging transition systems with high absolute initial h value",
        "0",
        Bounds("0", "100"));
    parser.add_option<int>(
        "w_high_absolute_max_f_value",
        "prefer merging transition systems with high absolute max f value",
        "0",
        Bounds("0", "100"));
    parser.add_option<int>(
        "w_high_absolute_max_g_value",
        "prefer merging transition systems with high absolute max g value",
        "0",
        Bounds("0", "100"));
    parser.add_option<int>(
        "w_high_absolute_max_h_value",
        "prefer merging transition systems with high absolute max h value",
        "0",
        Bounds("0", "100"));
    parser.add_option<int>(
        "w_high_average_h_value_improvement",
        "prefer merging transition systems with high average h value",
        "0",
        Bounds("0", "100"));
    parser.add_option<int>(
        "w_high_initial_h_value_sum",
        "prefer merging transition systems with large number of states",
        "0",
        Bounds("0", "100"));
    parser.add_option<int>(
        "w_high_average_h_value_sum",
        "prefer merging transition systems with large number of edges",
        "0",
        Bounds("0", "100"));
    parser.add_option<int>(
        "w_dfp",
        "merge according to DFP merge strategy",
        "0",
        Bounds("0", "100"));
    parser.add_option<int>(
        "w_goal_relevance",
        "prefer goal relevant transition systems",
        "0",
        Bounds("0", "100"));
    parser.add_option<int>(
        "w_num_variables",
        "prefer transition systems with many incorporated variables",
        "0",
        Bounds("0", "100"));
    parser.add_option<int>(
        "w_shrink_perfectly",
        "prefer merges which allow shrinking perfectly",
        "0",
        Bounds("0", "100"));
    parser.add_option<int>(
        "w_num_trans",
        "prefer transition systems with few transitions",
        "0",
        Bounds("0", "100"));
    parser.add_option<int>(
        "w_lr_opp",
        "prefer transition systems that allow for most label reductions",
        "0",
        Bounds("0", "100"));
    parser.add_option<int>(
        "w_more_lr_opp",
        "prefer transition systems that allow for most label reductions",
        "0",
        Bounds("0", "100"));
    parser.add_option<int>(
        "w_miasm",
        "prefer transition systems that allow for most unreachable and irrelevant pruning",
        "0",
        Bounds("0", "100"));
    parser.add_option<int>(
        "w_mutex",
        "prefer transition systems that have facts mutex to each other",
        "0",
        Bounds("0", "100"));

    Options opts = parser.parse();
    if (opts.get<int>("w_causally_connected_vars") == 0 &&
        opts.get<int>("w_bool_causally_connected_vars") == 0 &&
        opts.get<int>("w_nonadditive_vars") == 0 &&
        opts.get<int>("w_small_transitions_states_quotient") == 0 &&
        opts.get<int>("w_high_transitions_states_quotient") == 0 &&
        opts.get<int>("w_high_initial_h_value_improvement") == 0 &&
        opts.get<int>("w_high_absolute_initial_h_value") == 0 &&
        opts.get<int>("w_high_absolute_max_f_value") == 0 &&
        opts.get<int>("w_high_absolute_max_g_value") == 0 &&
        opts.get<int>("w_high_absolute_max_h_value") == 0 &&
        opts.get<int>("w_high_average_h_value_improvement") == 0 &&
        opts.get<int>("w_high_initial_h_value_sum") == 0 &&
        opts.get<int>("w_high_average_h_value_sum") == 0 &&
        opts.get<int>("w_dfp") == 0 &&
        opts.get<int>("w_goal_relevance") == 0 &&
        opts.get<int>("w_num_variables") == 0 &&
        opts.get<int>("w_shrink_perfectly") == 0 &&
        opts.get<int>("w_num_trans") == 0 &&
        opts.get<int>("w_lr_opp") == 0 &&
        opts.get<int>("w_more_lr_opp") == 0 &&
        opts.get<int>("w_miasm") == 0 &&
        opts.get<int>("w_mutex") == 0) {
        cerr << "you must specify at least one non-zero weight!" << endl;
        exit_with(EXIT_INPUT_ERROR);
    }

    if (parser.dry_run())
        return nullptr;
    else
        return make_shared<MergeDynamicWeighted>(opts);
}

static PluginShared<MergeStrategy> _plugin("merge_dynamic_weighted", _parse);
