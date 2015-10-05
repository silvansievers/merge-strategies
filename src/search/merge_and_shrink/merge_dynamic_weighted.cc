#include "merge_dynamic_weighted.h"

#include "transition_system.h"

#include "../causal_graph.h"
#include "../option_parser.h"
#include "../plugin.h"
#include "../task_proxy.h"


using namespace std;

// Helper methods to deal with transition systems

int compute_number_of_product_transitions(
    const TransitionSystem *ts1, const TransitionSystem *ts2) {
    // NOTE: this is copied from the merge constructor of TransitionSystem
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

AbstractFeature::AbstractFeature(bool merge_required)
    : merge_required(merge_required) {
}

CausalConnectionFeature::CausalConnectionFeature(const shared_ptr<AbstractTask> task_)
    : AbstractFeature(false),
      task(task_),
      causal_graph(TaskProxy(*task).get_causal_graph()) {
}

double CausalConnectionFeature::compute_value(const TransitionSystem *ts1,
                                              const TransitionSystem *ts2,
                                              const TransitionSystem *) {
    const vector<int> ts1_var_nos = ts1->get_incorporated_variables();
    vector<int> ts1_cg_neighbors;
    for (int var_no : ts1_var_nos) {
        const vector<int> &ts1_cg_successors = causal_graph.get_successors(var_no);
        ts1_cg_neighbors.insert(ts1_cg_neighbors.end(), ts1_cg_successors.begin(),
                                ts1_cg_successors.end());
        const vector<int> &ts1_cg_predecessors = causal_graph.get_predecessors(var_no);
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
    return static_cast<double>(edge_count) / max_possible_edges;
}

void CausalConnectionFeature::dump_precomputed_data() const {
    TaskProxy task_proxy(*task);
    cout << "causal graph:" << endl;
    for (VariableProxy var : task_proxy.get_variables()) {
        cout << "successors for var " << var.get_id() << ": "
             << causal_graph.get_successors(var.get_id()) << endl;
    }
}

NonAdditivityFeature::NonAdditivityFeature(const shared_ptr<AbstractTask> task)
    : AbstractFeature(false) {
    TaskProxy task_proxy(*task);
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
}

double NonAdditivityFeature::compute_value(const TransitionSystem *ts1,
                                           const TransitionSystem *ts2,
                                           const TransitionSystem *) {
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
    return static_cast<double>(not_additive_pair_count) / total_pair_count;
}

void NonAdditivityFeature::dump_precomputed_data() const {
    int num_variables = additive_var_pairs.size();
    for (int var_no1 = 0; var_no1 < num_variables; ++var_no1) {
        for (int var_no2 = var_no1 + 1; var_no2 < num_variables; ++var_no2) {
            cout << var_no1 << " and " << var_no2 << ": "
                 << (additive_var_pairs[var_no1][var_no2] ? "" : "not ")
                 << "additive" << endl;
        }
    }
}

TransStatesQuotFeature::TransStatesQuotFeature()
    : AbstractFeature(false) {
}

double TransStatesQuotFeature::compute_value(const TransitionSystem *ts1,
                                             const TransitionSystem *ts2,
                                             const TransitionSystem *) {
    double product_states = ts1->get_size() * ts2->get_size();
    double product_transitions = compute_number_of_product_transitions(ts1, ts2);
    return product_states / product_transitions;
}

InitHImprovementFeature::InitHImprovementFeature()
    : AbstractFeature(true) {
}

double InitHImprovementFeature::compute_value(const TransitionSystem *ts1,
                                              const TransitionSystem *ts2,
                                              const TransitionSystem *merge) {
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
        return 1;
    }
    return static_cast<double>(difference) / static_cast<double>(old_init_h);
}

AvgHImprovementFeature::AvgHImprovementFeature()
    : AbstractFeature(true) {
}

double AvgHImprovementFeature::compute_value(const TransitionSystem *ts1,
                                             const TransitionSystem *ts2,
                                             const TransitionSystem *merge) {
    assert(merge);
    double new_average_h = compute_average_h_value(merge);
    double old_average_h = max(compute_average_h_value(ts1),
                               compute_average_h_value(ts2));
    double difference = new_average_h - old_average_h;
    if (!difference) {
        return 0;
    }
    if (!old_average_h) {
        return 1;
    }
    return static_cast<double>(difference) / static_cast<double>(old_average_h);
}

InitHSumFeature::InitHSumFeature()
    : AbstractFeature(false) {
}

double InitHSumFeature::compute_value(const TransitionSystem *ts1,
                                      const TransitionSystem *ts2,
                                      const TransitionSystem *) {
    int init_h_sum = ts1->get_init_state_goal_distance() +
        ts2->get_init_state_goal_distance();
    return init_h_sum;
}

AvgHSumFeature::AvgHSumFeature()
    : AbstractFeature(false) {
}

double AvgHSumFeature::compute_value(const TransitionSystem *ts1,
                                     const TransitionSystem *ts2,
                                     const TransitionSystem *) {
    double average_h_sum = compute_average_h_value(ts1) +
        compute_average_h_value(ts2);
    return average_h_sum;
}

Features::Features(std::vector<int> &&weights_,
                   const shared_ptr<AbstractTask> task,
                   bool debug)
    : weights(move(weights_)),
      task(task),
      debug(debug),
      features(weights.size(), 0) {
    clear();
    if (weights[0]) {
        features[0] = new CausalConnectionFeature(task);
    }
    if (weights[1]) {
        features[1] = new NonAdditivityFeature(task);
    }
    if (weights[2]) {
        features[2] = new TransStatesQuotFeature();
    }
    if (weights[3]) {
        features[3] = new InitHImprovementFeature();
    }
    if (weights[4]) {
        features[4] = new AvgHImprovementFeature();
    }
    if (weights[5]) {
        features[5] = new InitHSumFeature();
    }
    if (weights[6]) {
        features[6] = new AvgHSumFeature();
    }
    if (debug) {
        for (AbstractFeature *feature : features) {
            if (feature) {
                feature->dump_precomputed_data();
            }
        }
    }
}

Features::~Features() {
    for (AbstractFeature *feature : features) {
        delete feature;
    }
}

void Features::update_min_max(int feature_no, double value) {
    if (value > max_values[feature_no]) {
        max_values[feature_no] = value;
    }
    if (value < min_values[feature_no]) {
        min_values[feature_no] = value;
    }
}

double Features::normalize_value(int feature_no, double value) const {
    double min = min_values[feature_no];
    double max = max_values[feature_no];
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
                                              TransitionSystem *ts2) {
    TransitionSystem *merge = 0;
    vector<double> values;
    values.reserve(features.size());
    for (size_t feature_no = 0; feature_no < features.size(); ++feature_no) {
        AbstractFeature *feature = features[feature_no];
        if (feature) {
            if (feature->requires_merge() && !merge) {
                merge = new TransitionSystem(TaskProxy(*task),
                                             ts1->get_labels(),
                                             ts1, ts2, false);
            }
            double value = feature->compute_value(ts1, ts2, merge);
            update_min_max(feature_no, value);
            values.push_back(value);
        } else {
            // dummy value for correct indices
            values.push_back(-1);
        }
    }
    unnormalized_values[make_pair(ts1, ts2)] = values;
    delete merge;
}

double Features::compute_weighted_normalized_sum(
    TransitionSystem *ts1, TransitionSystem *ts2) const {
    const std::vector<double> &values = unnormalized_values.at(make_pair(ts1, ts2));
    double weighted_sum = 0;
    if (debug) {
        cout << "computing weighted normalized sum for "
             << ts1->tag() << ts2->tag() << endl;
    }
    for (size_t feature_no = 0; feature_no < features.size(); ++feature_no) {
        if (weights[feature_no]) {
            double normalized_value = normalize_value(feature_no, values[feature_no]);
            if (debug) {
                cout << "normalized value for feature number " << feature_no
                     << ": " << normalized_value << endl;
            }
            weighted_sum += weights[feature_no] * normalized_value;
        }
    }
    if (debug) {
        cout << "weighted normalized sum: " << weighted_sum << endl;
    }
    return weighted_sum;
}

void Features::clear() {
    min_values.assign(features.size(), INF);
    max_values.assign(features.size(), -1);
    unnormalized_values.clear();
}

MergeDynamicWeighted::MergeDynamicWeighted(const Options opts)
    : MergeStrategy(),
      debug(opts.get<bool>("debug")),
      w_causally_connected_vars(opts.get<int>("w_causally_connected_vars")),
      w_nonadditive_vars(opts.get<int>("w_nonadditive_vars")),
      w_small_transitions_states_quotient(opts.get<int>("w_small_transitions_states_quotient")),
      w_high_initial_h_value_improvement(opts.get<int>("w_high_initial_h_value_improvement")),
      w_high_average_h_value_improvement(opts.get<int>("w_high_average_h_value_improvement")),
      w_high_initial_h_value_sum(opts.get<int>("w_high_initial_h_value_sum")),
      w_high_average_h_value_sum(opts.get<int>("w_high_average_h_value_sum")),
      features(0) {
}

MergeDynamicWeighted::~MergeDynamicWeighted() {
    delete features;
}

void MergeDynamicWeighted::dump_strategy_specific_options() const {
    cout << "w_causally_connected_vars: " << w_causally_connected_vars << endl;
    cout << "w_nonadditive_vars: " << w_nonadditive_vars << endl;
    cout << "w_small_transitions_states_quotient: " << w_small_transitions_states_quotient << endl;
    cout << "w_high_initial_h_value_improvement: " << w_high_initial_h_value_improvement << endl;
    cout << "w_high_average_h_value_improvement: " << w_high_average_h_value_improvement << endl;
    cout << "w_high_initial_h_value_sum: " << w_high_initial_h_value_sum << endl;
    cout << "w_high_average_h_value_sum: " << w_high_average_h_value_sum << endl;
}

void MergeDynamicWeighted::initialize(const shared_ptr<AbstractTask> task) {
    MergeStrategy::initialize(task);
    TaskProxy task_proxy(*task);
    int num_variables = task_proxy.get_variables().size();
    var_no_to_ts_index.reserve(num_variables);
    for (VariableProxy var : task_proxy.get_variables()) {
        var_no_to_ts_index.push_back(var.get_id());
    }
    merge_order.reserve(num_variables * 2 - 1);
    vector<int> feature_weights = {
        w_causally_connected_vars,
        w_nonadditive_vars,
        w_small_transitions_states_quotient,
        w_high_initial_h_value_improvement,
        w_high_average_h_value_improvement,
        w_high_initial_h_value_sum,
        w_high_average_h_value_sum
    };
    features = new Features(move(feature_weights), task, debug);
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
        // Go through all transitition systems and compute unnormalized feature values.
        for (size_t i = 0; i < all_transition_systems.size(); ++i) {
            TransitionSystem *ts1 = all_transition_systems[i];
            if (ts1) {
                for (size_t j = i + 1; j < all_transition_systems.size(); ++j) {
                    TransitionSystem *ts2 = all_transition_systems[j];
                    if (ts2) {
                        features->precompute_unnormalized_values(ts1, ts2);
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
                            features->compute_weighted_normalized_sum(
                                all_transition_systems[i], all_transition_systems[j]);
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

    Options opts = parser.parse();
    if (opts.get<int>("w_causally_connected_vars") == 0 &&
        opts.get<int>("w_nonadditive_vars") == 0 &&
        opts.get<int>("w_small_transitions_states_quotient") == 0 &&
        opts.get<int>("w_high_initial_h_value_improvement") == 0 &&
        opts.get<int>("w_high_average_h_value_improvement") == 0 &&
        opts.get<int>("w_high_initial_h_value_sum") == 0 &&
        opts.get<int>("w_high_average_h_value_sum") == 0) {
        cerr << "you must specify at least one non-zero weight!" << endl;
        exit_with(EXIT_INPUT_ERROR);
    }

    if (parser.dry_run())
        return nullptr;
    else
        return make_shared<MergeDynamicWeighted>(opts);
}

static PluginShared<MergeStrategy> _plugin("merge_dynamic_weighted", _parse);
