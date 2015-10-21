#include "merge_dynamic_weighted.h"

#include "labels.h"
#include "shrink_bisimulation.h"
#include "transition_system.h"

#include "../causal_graph.h"
#include "../option_parser.h"
#include "../plugin.h"
#include "../task_proxy.h"


using namespace std;

const int MINUSINF = numeric_limits<int>::min();

// Helper methods to deal with transition systems

int compute_number_of_product_transitions(
    const TransitionSystem *ts1, const TransitionSystem *ts2) {
    // TODO: this is copied from the merge constructor of TransitionSystem
    /*
      Note that this computes the number of tranistions in the product
      without considering possible shrinking due to unreachable or
      irrelevant states, which hence may reduce the actual number of
      transitions in the product.
    */
    int number_of_transitions = 0;
    for (TSConstIterator group1_it = ts1->begin();
         group1_it != ts1->end(); ++group1_it) {
        // Distribute the labels of this group among the "buckets"
        // corresponding to the groups of ts2.
        unordered_map<int, vector<int> > buckets;
        for (LabelConstIter label_it = group1_it.begin();
             label_it != group1_it.end(); ++label_it) {
            int label_no = *label_it;
            int group2_id = ts2->get_group_id_for_label(label_no);
            buckets[group2_id].push_back(label_no);
        }
        // Now buckets contains all equivalence classes that are
        // refinements of group1.

        // Now create the new groups together with their transitions.
        const vector<Transition> &transitions1 = group1_it.get_transitions();
        for (const auto &bucket : buckets) {
            const vector<Transition> &transitions2 =
                ts2->get_transitions_for_group_id(bucket.first);
            int new_transitions_for_new_group = transitions1.size() * transitions2.size();
            number_of_transitions += new_transitions_for_new_group;
        }
    }
    return number_of_transitions;
}

double compute_average_h_value(const TransitionSystem *ts) {
    int num_states = ts->get_size();
    int sum_distances = 0;
    for (int state = 0; state < num_states; ++state) {
        sum_distances += ts->get_goal_distance(state);
    }
    if (num_states == 0) {
        // For unsolvable transition systems
        return INF;
    }
    return static_cast<double>(sum_distances) / static_cast<double>(num_states);
}

void compute_label_ranks(const TransitionSystem *ts,
                         vector<int> &label_ranks) {
    // TODO: copied from MergeDFP
    int num_labels = ts->get_num_labels();
    // Irrelevant (and inactive, i.e. reduced) labels have a dummy rank of -1
    label_ranks.resize(num_labels, -1);

    for (TSConstIterator group_it = ts->begin();
         group_it != ts->end(); ++group_it) {
        // Relevant labels with no transitions have a rank of infinity.
        int label_rank = INF;
        const vector<Transition> &transitions = group_it.get_transitions();
        bool group_relevant = false;
        if (static_cast<int>(transitions.size()) == ts->get_size()) {
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
                label_rank = min(label_rank, ts->get_goal_distance(t.target));
            }
        }
        for (LabelConstIter label_it = group_it.begin();
             label_it != group_it.end(); ++label_it) {
            int label_no = *label_it;
            label_ranks[label_no] = label_rank;
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

double Feature::compute_unnormalized_value(TransitionSystem *ts1,
                                           TransitionSystem *ts2,
                                           TransitionSystem *merge) {
    if (weight) {
        // TODO: get rid of this? we currently also perform this check
        // outside of this method before calling it
        return compute_value(ts1, ts2, merge);
    }
    return 0;
}

void Feature::dump() const {
    cout << name << ": " << weight << endl;
}

CausalConnectionFeature::CausalConnectionFeature(int id, int weight)
    : Feature(id, weight, "causally connected variables", false, false),
      causal_graph(0) {
}

CausalConnectionFeature::~CausalConnectionFeature() {
    delete causal_graph;
}

double CausalConnectionFeature::compute_value(TransitionSystem *ts1,
                                              TransitionSystem *ts2,
                                              TransitionSystem *) {
    // return value in [0,infinity)
    const vector<int> ts1_var_nos = ts1->get_incorporated_variables();
    vector<int> ts1_cg_neighbors;
    for (int var_no : ts1_var_nos) {
        const vector<int> &ts1_cg_successors = causal_graph->get_successors(var_no);
        ts1_cg_neighbors.insert(ts1_cg_neighbors.end(), ts1_cg_successors.begin(),
                                ts1_cg_successors.end());
        const vector<int> &ts1_cg_predecessors = causal_graph->get_predecessors(var_no);
        ts1_cg_neighbors.insert(ts1_cg_neighbors.end(), ts1_cg_predecessors.begin(),
                                ts1_cg_predecessors.end());
    }
    // NOTE: we want to count directed edges, to take into account if there
    // is a connection between variables in both directions or not.
//    sort(ts1_cg_neighbors.begin(), ts1_cg_neighbors.end());
//    ts1_cg_neighbors.erase(unique(ts1_cg_neighbors.begin(), ts1_cg_neighbors.end()),
//                           ts1_cg_neighbors.end());

    const vector<int> ts2_var_nos = ts2->get_incorporated_variables();
    int edge_count = 0;
    for (int ts1_cg_neighbor : ts1_cg_neighbors) {
        for (int ts2_var_no : ts2_var_nos) {
            if (ts2_var_no == ts1_cg_neighbor) {
                ++edge_count;
            }
        }
    }
    double max_possible_edges = ts1_var_nos.size() * ts2_var_nos.size() * 2;
    assert(max_possible_edges);
    return static_cast<double>(edge_count) / max_possible_edges;
}

void CausalConnectionFeature::initialize(const TaskProxy &task_proxy, bool dump) {
    // TODO: avoid recreating a new causal graph (cannot assign a const
    // reference in this method)
    causal_graph = new CausalGraph(task_proxy);
    if (dump) {
        cout << "causal graph:" << endl;
        for (VariableProxy var : task_proxy.get_variables()) {
            cout << "successors for var " << var.get_id() << ": "
                 << causal_graph->get_successors(var.get_id()) << endl;
        }
    }
}

NonAdditivityFeature::NonAdditivityFeature(int id, int weight)
    : Feature(id, weight, "non additive variables", false, false) {
}

double NonAdditivityFeature::compute_value(TransitionSystem *ts1,
                                           TransitionSystem *ts2,
                                           TransitionSystem *) {
    // return value in [0,infinity)
    const vector<int> ts1_var_nos = ts1->get_incorporated_variables();
    const vector<int> ts2_var_nos = ts2->get_incorporated_variables();
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

double SmallTransStatesQuotFeature::compute_value(TransitionSystem *ts1,
                                                  TransitionSystem *ts2,
                                                  TransitionSystem *) {
    // return value in [0,infinity)
    double product_states = ts1->get_size() * ts2->get_size();
    double product_transitions = compute_number_of_product_transitions(ts1, ts2);
    if (!product_states) {
        return INF;
    }
    return product_transitions / product_states;
}

HighTransStatesQuotFeature::HighTransStatesQuotFeature(int id, int weight)
    : Feature(id, weight, "high transitions per states quotient", false, false) {
}

double HighTransStatesQuotFeature::compute_value(TransitionSystem *ts1,
                                                 TransitionSystem *ts2,
                                                 TransitionSystem *) {
    // return value in [0,infinity)
    double product_states = ts1->get_size() * ts2->get_size();
    double product_transitions = compute_number_of_product_transitions(ts1, ts2);
    if (!product_states) {
        return INF;
    }
    return product_transitions / product_states;
}

InitHImprovementFeature::InitHImprovementFeature(int id, int weight)
    : Feature(id, weight, "initial h value improvement", true, false) {
}

double InitHImprovementFeature::compute_value(TransitionSystem *ts1,
                                              TransitionSystem *ts2,
                                              TransitionSystem *merge) {
    // return value in [0,infinity)
    assert(merge);
    int new_init_h;
    if (merge->is_solvable()) {
        new_init_h = merge->get_init_state_goal_distance();
    } else {
        // initial state has been pruned
        new_init_h = INF;
    }
    int old_init_h = max(ts1->get_init_state_goal_distance(),
                         ts2->get_init_state_goal_distance());
    int difference = new_init_h - old_init_h;
    if (!difference) {
        return 0;
    }
    if (!old_init_h) {
        return INF;
    }
    return static_cast<double>(difference) / static_cast<double>(old_init_h);
}

AvgHImprovementFeature::AvgHImprovementFeature(int id, int weight)
    : Feature(id, weight, "average h value improvement", true, false) {
}

double AvgHImprovementFeature::compute_value(TransitionSystem *ts1,
                                             TransitionSystem *ts2,
                                             TransitionSystem *merge) {
    // return value in [0,infinity)
    assert(merge);
    double new_average_h = compute_average_h_value(merge);
    double old_average_h = max(compute_average_h_value(ts1),
                               compute_average_h_value(ts2));
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

double InitHSumFeature::compute_value(TransitionSystem *ts1,
                                      TransitionSystem *ts2,
                                      TransitionSystem *) {
    // return value in [0,infinity)
    int init_h_sum = ts1->get_init_state_goal_distance() +
        ts2->get_init_state_goal_distance();
    return init_h_sum;
}

AvgHSumFeature::AvgHSumFeature(int id, int weight)
    : Feature(id, weight, "average h value sum", false, false) {
}

double AvgHSumFeature::compute_value(TransitionSystem *ts1,
                                     TransitionSystem *ts2,
                                     TransitionSystem *) {
    // return value in [0,infinity)
    double average_h_sum = compute_average_h_value(ts1) +
        compute_average_h_value(ts2);
    return average_h_sum;
}

DFPFeature::DFPFeature(int id, int weight)
    : Feature(id, weight, "dfp", false, true) {
}

void DFPFeature::clear() {
    ts_to_label_ranks.clear();
}

double DFPFeature::compute_value(TransitionSystem *ts1,
                                 TransitionSystem *ts2,
                                 TransitionSystem *) {
    // return value in [0,infinity)
    int pair_weight = INF;
    if (ts1->is_goal_relevant() || ts2->is_goal_relevant()) {
        vector<int> &label_ranks1 = ts_to_label_ranks[ts1];
        if (label_ranks1.empty()) {
            compute_label_ranks(ts1, label_ranks1);
        }
        vector<int> &label_ranks2 = ts_to_label_ranks[ts2];
        if (label_ranks2.empty()) {
            compute_label_ranks(ts2, label_ranks2);
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

double GoalRelevanceFeature::compute_value(TransitionSystem *ts1,
                                           TransitionSystem *ts2,
                                           TransitionSystem *) {
    // return value in [0,2]
    int pair_weight = 0;
    if (ts1->is_goal_relevant()) {
        ++pair_weight;
    }
    if (ts2->is_goal_relevant()) {
        ++pair_weight;
    }
    return pair_weight;
}

NumVariablesFeature::NumVariablesFeature(int id, int weight)
    : Feature(id, weight, "high number of incorporated variables", false, false) {
}

double NumVariablesFeature::compute_value(TransitionSystem *ts1,
                                          TransitionSystem *ts2,
                                          TransitionSystem *) {
    // return value in [2,num_variables-1]
    return ts1->get_incorporated_variables().size() +
        ts2->get_incorporated_variables().size();
}

NumTransitionsFeature::NumTransitionsFeature(int id, int weight)
    : Feature(id, weight, "small number of transitions", false, true) {
}

double NumTransitionsFeature::compute_value(TransitionSystem *ts1,
                                            TransitionSystem *ts2,
                                            TransitionSystem *) {
    // return value in [0,infinity[
    return compute_number_of_product_transitions(ts1, ts2);
}

LROpportunitiesFeatures::LROpportunitiesFeatures(int id, int weight)
    : Feature(id, weight, "high number of reducible labels", false, false) {
}

void LROpportunitiesFeatures::clear() {
    ts_pair_to_combinable_label_count.clear();
}

void LROpportunitiesFeatures::precompute_data(
    const vector<TransitionSystem *> &all_transition_systems) {
    // Precompute the set of irrelevant labels for every transition system
    unordered_map<TransitionSystem *, vector<bool> > ts_to_irrelevant_labels;
    for (TransitionSystem *ts : all_transition_systems) {
        if (ts) {
            vector<bool> irrelevant_labels(ts->get_num_labels(), false);
            for (TSConstIterator group_it = ts->begin();
                 group_it != ts->end(); ++group_it) {
                const vector<Transition> &transitions = group_it.get_transitions();
                bool group_relevant = false;
                if (static_cast<int>(transitions.size()) == ts->get_size()) {
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
            ts_to_irrelevant_labels[ts] = irrelevant_labels;
        }
    }

    // Compute the number of labels that are irrelevant in all other transition
    // systems than then current considered pair.
    const shared_ptr<Labels> labels = all_transition_systems.back()->get_labels();
    int num_labels = labels->get_size();
    for (size_t i = 0; i < all_transition_systems.size(); ++i) {
        TransitionSystem *ts1 = all_transition_systems[i];
        if (ts1) {
            for (size_t j = i + 1; j < all_transition_systems.size(); ++j) {
                TransitionSystem *ts2 = all_transition_systems[j];
                if (ts2) {
                    int count_combinable_labels = 0;
                    for (int label_no = 0; label_no < num_labels; ++label_no) {
                        bool label_irrelevant_in_all_other_ts = true;
                        for (size_t k = 0; k < all_transition_systems.size(); ++k) {
                            if (k == i || k == j) {
                                continue;
                            }
                            TransitionSystem *ts3 = all_transition_systems[k];
                            if (ts3) {
                                if (!ts_to_irrelevant_labels[ts3][label_no]) {
                                    label_irrelevant_in_all_other_ts = false;
                                    break;
                                }
                            }

                        }
                        if (label_irrelevant_in_all_other_ts) {
                            ++count_combinable_labels;
                        }
                    }
                    ts_pair_to_combinable_label_count[make_pair(ts1, ts2)] =
                        count_combinable_labels;
                }
            }
        }
    }
}

double LROpportunitiesFeatures::compute_value(TransitionSystem *ts1,
                                              TransitionSystem *ts2,
                                              TransitionSystem *) {
    // return value in [0,infinity[
    int combinable_labels = ts_pair_to_combinable_label_count[make_pair(ts1, ts2)];
    if (combinable_labels >= 2) {
        // need at least two labels to profit from label reduction
        return combinable_labels;
    }
    return 0;
}

// ========================= FEATURES ====================================

Features::Features(const Options opts)
    : debug(opts.get<bool>("debug")),
      merge_required(false) {
    int id = 0;
    features.push_back(new CausalConnectionFeature(
                           id++, opts.get<int>("w_causally_connected_vars")));
    features.push_back(new NonAdditivityFeature(
                           id++, opts.get<int>("w_nonadditive_vars")));
    features.push_back(new SmallTransStatesQuotFeature(
                           id++, opts.get<int>("w_small_transitions_states_quotient")));
    features.push_back(new HighTransStatesQuotFeature(
                           id++, opts.get<int>("w_high_transitions_states_quotient")));
    features.push_back(new InitHImprovementFeature(
                           id++, opts.get<int>("w_high_initial_h_value_improvement")));
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
    features.push_back(new NumTransitionsFeature(
                           id++, opts.get<int>("w_num_trans")));
    features.push_back(new LROpportunitiesFeatures(
                           id++, opts.get<int>("w_lr_opp")));
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
    const vector<TransitionSystem *> &all_transition_sytems) {
    for (Feature *feature : features) {
        if (feature->get_weight()) {
            feature->precompute_data(all_transition_sytems);
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

void Features::precompute_unnormalized_values(TransitionSystem *ts1,
                                              TransitionSystem *ts2,
                                              TransitionSystem *merge) {
    vector<double> values;
    values.reserve(features.size());
    for (Feature *feature : features) {
        if (feature->get_weight()) {
            double value = feature->compute_unnormalized_value(ts1, ts2, merge);
            update_min_max(feature->get_id(), value);
            values.push_back(value);
        } else {
            // dummy value for correct indices
            // (will be multiplied with weight 0 anyway)
            values.push_back(0);
        }
    }
    unnormalized_values[make_pair(ts1, ts2)] = values;
}

double Features::compute_weighted_normalized_sum(
    TransitionSystem *ts1, TransitionSystem *ts2) const {
    const std::vector<double> &values = unnormalized_values.at(make_pair(ts1, ts2));
    double weighted_sum = 0;
    if (debug) {
        cout << "computing weighted normalized sum for "
             << ts1->tag() << ts2->tag() << endl;
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
      max_states(opts.get<int>("max_states")) {
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
    int num_variables = task_proxy->get_variables().size();
    var_no_to_ts_index.reserve(num_variables);
    for (VariableProxy var : task_proxy->get_variables()) {
        var_no_to_ts_index.push_back(var.get_id());
    }
    merge_order.reserve(num_variables * 2 - 1);
    features->initialize(*task_proxy);
}

pair<int, int> MergeDynamicWeighted::get_next(const vector<TransitionSystem *> &all_transition_systems) {
    int ts_index1 = -1;
    int ts_index2 = -1;

    if (remaining_merges == 1) {
        // Return the first pair
        ts_index1 = 0;
        while (!all_transition_systems[ts_index1]) {
            ++ts_index1;
        }
        assert(in_bounds(ts_index1, all_transition_systems));
        ts_index2 = ts_index1 + 1;
        while (!all_transition_systems[ts_index2]) {
            ++ts_index2;
        }
        assert(in_bounds(ts_index2, all_transition_systems));
    } else {
        features->precompute_data(all_transition_systems);
        // Go through all transitition systems and compute unnormalized feature values.
        vector<TransitionSystem *> tmp_transition_systems(all_transition_systems);
        for (size_t i = 0; i < all_transition_systems.size(); ++i) {
            TransitionSystem *ts1 = all_transition_systems[i];
            if (ts1) {
                for (size_t j = i + 1; j < all_transition_systems.size(); ++j) {
                    TransitionSystem *ts2 = all_transition_systems[j];
                    if (ts2) {
                        TransitionSystem *merge = 0;
                        if (features->require_merge()) {
                            const shared_ptr<Labels> labels = ts1->get_labels();
                            shared_ptr<Labels> labels_copy = make_shared<Labels>(*labels.get());

                            TransitionSystem *ts1_copy = new TransitionSystem(*ts1, labels_copy);
                            TransitionSystem *ts2_copy = new TransitionSystem(*ts2, labels_copy);
                            tmp_transition_systems[i] = ts1_copy;
                            tmp_transition_systems[j] = ts2_copy;

                            if (labels_copy->reduce_before_shrinking()) {
                                labels_copy->reduce(make_pair(i, j),
                                                    tmp_transition_systems,
                                                    true);
                            }

                            Options options;
                            options.set<int>("max_states", max_states);
                            options.set<int>("max_states_before_merge", max_states);
                            options.set<int>("threshold", 1);
                            options.set<bool>("greedy", false);
                            options.set<int>("at_limit", 0);
                            ShrinkBisimulation shrink_bisim(options);
                            shrink_bisim.shrink(*ts1_copy, *ts2_copy, true);
                            merge = new TransitionSystem(*task_proxy,
                                                         labels_copy,
                                                         ts1_copy, ts2_copy, true);
                        }
                        features->precompute_unnormalized_values(ts1, ts2, merge);
                        if (features->require_merge()) {
                            // delete and reset
                            delete merge;
                            delete tmp_transition_systems[i];
                            delete tmp_transition_systems[j];
                            tmp_transition_systems[i] = ts1;
                            tmp_transition_systems[j] = ts2;
                        }
                    }
                }
            }
        }

        // Go through all transition systems again and normalize feature values.
        int max_weight = -1;
        for (size_t i = 0; i < all_transition_systems.size(); ++i) {
            TransitionSystem *ts1 = all_transition_systems[i];
            if (ts1) {
                for (size_t j = i + 1; j < all_transition_systems.size(); ++j) {
                    TransitionSystem *ts2 = all_transition_systems[j];
                    if (ts2) {
                        int pair_weight =
                            features->compute_weighted_normalized_sum(ts1, ts2);
                        if (pair_weight > max_weight) {
                            max_weight = pair_weight;
                            ts_index1 = i;
                            ts_index2 = j;
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

    assert(ts_index1 != -1);
    assert(ts_index2 != -1);
    int new_ts_index = all_transition_systems.size();
    TransitionSystem *ts1 = all_transition_systems[ts_index1];
    TransitionSystem *ts2 = all_transition_systems[ts_index2];
    for (int var_no : ts1->get_incorporated_variables()) {
        var_no_to_ts_index[var_no] = new_ts_index;
    }
    for (int var_no : ts2->get_incorporated_variables()) {
        var_no_to_ts_index[var_no] = new_ts_index;
    }
    cout << "Next pair of indices: (" << ts_index1 << ", " << ts_index2 << ")" << endl;
    --remaining_merges;
    merge_order.push_back(make_pair(ts_index1, ts_index2));
    if (!remaining_merges) {
        cout << "merge order: ";
        for (auto merge : merge_order) {
            cout << "(" << merge.first << ", " << merge.second << "), ";
        }
        cout << endl;
    }
    return make_pair(ts_index1, ts_index2);
}

std::string MergeDynamicWeighted::name() const {
    return "dynamic merging";
}

static shared_ptr<MergeStrategy>_parse(OptionParser &parser) {
    parser.add_option<bool>(
        "debug", "debug", "false");
    parser.add_option<int>("max_states", "shrink strategy option", "50000");
    parser.add_option<int>(
        "w_causally_connected_vars",
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
        "prefer merging transition systems with high initial h value",
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
        "w_num_trans",
        "prefer transition systems with few transitions",
        "0",
        Bounds("0", "100"));
    parser.add_option<int>(
        "w_lr_opp",
        "prefer transition systems that allow for most label reductions",
        "0",
        Bounds("0", "100"));

    Options opts = parser.parse();
    if (opts.get<int>("w_causally_connected_vars") == 0 &&
        opts.get<int>("w_nonadditive_vars") == 0 &&
        opts.get<int>("w_small_transitions_states_quotient") == 0 &&
        opts.get<int>("w_high_transitions_states_quotient") == 0 &&
        opts.get<int>("w_high_initial_h_value_improvement") == 0 &&
        opts.get<int>("w_high_average_h_value_improvement") == 0 &&
        opts.get<int>("w_high_initial_h_value_sum") == 0 &&
        opts.get<int>("w_high_average_h_value_sum") == 0 &&
        opts.get<int>("w_dfp") == 0 &&
        opts.get<int>("w_goal_relevance") == 0 &&
        opts.get<int>("w_num_variables") == 0 &&
        opts.get<int>("w_num_trans") == 0 &&
        opts.get<int>("w_lr_opp") == 0) {
        cerr << "you must specify at least one non-zero weight!" << endl;
        exit_with(EXIT_INPUT_ERROR);
    }

    if (parser.dry_run())
        return nullptr;
    else
        return make_shared<MergeDynamicWeighted>(opts);
}

static PluginShared<MergeStrategy> _plugin("merge_dynamic_weighted", _parse);
