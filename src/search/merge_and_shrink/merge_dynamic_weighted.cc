#include "merge_dynamic_weighted.h"

#include "transition_system.h"

#include "../causal_graph.h"
#include "../option_parser.h"
#include "../plugin.h"
#include "../task_proxy.h"


using namespace std;

MergeDynamicWeighted::MergeDynamicWeighted(const Options opts)
    : MergeStrategy(),
      w_prefer_causally_connected_vars(opts.get<int>("w_prefer_causally_connected_vars")),
      w_avoid_additive_vars(opts.get<int>("w_avoid_additive_vars")),
      w_high_initial_h_value(opts.get<int>("w_high_initial_h_value")),
      w_high_average_h_value(opts.get<int>("w_high_average_h_value")),
      w_prefer_ts_large_num_states(opts.get<int>("w_prefer_ts_large_num_states")),
      w_prefer_ts_large_num_edges(opts.get<int>("w_prefer_ts_large_num_edges")),
      causal_graph(0) {
}

MergeDynamicWeighted::~MergeDynamicWeighted() {
    delete causal_graph;
}

void MergeDynamicWeighted::dump_strategy_specific_options() const {
    cout << "w_prefer_causally_connected_vars: " << w_prefer_causally_connected_vars << endl;
    cout << "w_avoid_additive_vars: " << w_avoid_additive_vars << endl;
    cout << "w_high_initial_h_value: " << w_high_initial_h_value << endl;
    cout << "w_high_average_h_value: " << w_high_average_h_value << endl;
    cout << "w_prefer_ts_large_num_states: " << w_prefer_ts_large_num_states << endl;
    cout << "w_prefer_ts_large_num_edges: " << w_prefer_ts_large_num_edges << endl;
}

void MergeDynamicWeighted::initialize(const shared_ptr<AbstractTask> task_) {
    MergeStrategy::initialize(task_);
    task = task_;
    TaskProxy task_proxy(*task);
    int num_variables = task_proxy.get_variables().size();
    var_no_to_ts_index.reserve(num_variables);
    for (VariableProxy var : task_proxy.get_variables()) {
        var_no_to_ts_index.push_back(var.get_id());
    }
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
    if (w_prefer_causally_connected_vars) {
        // TODO: do not recreate causal graph. This is a solution for
        // circumventing of assigning a const reference to a non-const member.
        causal_graph = new CausalGraph(task_proxy.get_causal_graph());
    }
}

int MergeDynamicWeighted::compute_pair_weight(
    TransitionSystem *ts1, TransitionSystem *ts2) const {
    int weight = -1;
    if (w_prefer_causally_connected_vars) {
        const vector<int> ts1_var_nos = ts1->get_incorporated_variables();
        vector<int> ts1_cg_neighbors;
        for (int var_no : ts1_var_nos) {
            const vector<int> &ts1_cg_successors = causal_graph->get_successors(var_no);
            ts1_cg_neighbors.insert(ts1_cg_neighbors.end(), ts1_cg_successors.begin(), ts1_cg_successors.end());
            const vector<int> &ts1_cg_predecessors = causal_graph->get_predecessors(var_no);
            ts1_cg_neighbors.insert(ts1_cg_neighbors.end(), ts1_cg_predecessors.begin(), ts1_cg_predecessors.end());
        }
        sort(ts1_cg_neighbors.begin(), ts1_cg_neighbors.end());
        ts1_cg_neighbors.erase(unique(ts1_cg_neighbors.begin(), ts1_cg_neighbors.end()),
                               ts1_cg_neighbors.end());

        const vector<int> ts2_var_nos = ts2->get_incorporated_variables();
        bool causally_connected = false;
        for (int ts1_cg_neighbor : ts1_cg_neighbors) {
            for (int ts2_var_no : ts2_var_nos) {
                if (ts2_var_no == ts1_cg_neighbor) {
                    causally_connected = true;
                    break;
                }
            }
        }
        if (causally_connected) {
            weight += w_prefer_causally_connected_vars;
        }
    }
    return weight;
}

pair<int, int> MergeDynamicWeighted::get_next(const vector<TransitionSystem *> &all_transition_systems) {
    int ts_index1 = -1;
    int ts_index2 = -1;
    int max_weight = -1;
    for (size_t i = 0; i < all_transition_systems.size(); ++i) {
        TransitionSystem *ts1 = all_transition_systems[i];
        if (ts1) {
            for (size_t j = i + 1; j < all_transition_systems.size(); ++j) {
                TransitionSystem *ts2 = all_transition_systems[j];
                if (ts2) {
                    int pair_weight = compute_pair_weight(all_transition_systems[i], all_transition_systems[j]);
                    if (pair_weight > max_weight) {
                        max_weight = pair_weight;
                        ts_index1 = i;
                        ts_index2 = j;
                    }
                }
            }
        }
    }

    if (max_weight == -1) {
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
    }

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
    return make_pair(ts_index1, ts_index2);
}

std::string MergeDynamicWeighted::name() const {
    return "dynamic merging";
}

static shared_ptr<MergeStrategy>_parse(OptionParser &parser) {
    parser.add_option<int>(
        "w_prefer_causally_connected_vars",
        "prefer merging variables that are causally connected ",
        "0",
        Bounds("0", "100"));
    parser.add_option<int>(
        "w_avoid_additive_vars",
        "avoid merging additive variables",
        "0",
        Bounds("0", "100"));
    parser.add_option<int>(
        "w_high_initial_h_value",
        "prefer merging transition systems with high initial h value",
        "0",
        Bounds("0", "100"));
    parser.add_option<int>(
        "w_high_average_h_value",
        "prefer merging transition systems with high average h value",
        "0",
        Bounds("0", "100"));
    parser.add_option<int>(
        "w_prefer_ts_large_num_states",
        "prefer merging transition systems with large number of states",
        "0",
        Bounds("0", "100"));
    parser.add_option<int>(
        "w_prefer_ts_large_num_edges",
        "prefer merging transition systems with large number of edges",
        "0",
        Bounds("0", "100"));

    Options opts = parser.parse();
    if (opts.get<int>("w_prefer_causally_connected_vars") == 0 &&
        opts.get<int>("w_avoid_additive_vars") == 0 &&
        opts.get<int>("w_high_initial_h_value") == 0 &&
        opts.get<int>("w_high_average_h_value") == 0 &&
        opts.get<int>("w_prefer_ts_large_num_states") == 0 &&
        opts.get<int>("w_prefer_ts_large_num_edges") == 0) {
        cerr << "you must specify at least one non-zero weight!" << endl;
        exit_with(EXIT_INPUT_ERROR);
    }

    if (parser.dry_run())
        return nullptr;
    else
        return make_shared<MergeDynamicWeighted>(opts);
}

static PluginShared<MergeStrategy> _plugin("merge_dynamic_weighted", _parse);
