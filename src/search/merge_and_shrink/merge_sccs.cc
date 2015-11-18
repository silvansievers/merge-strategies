#include "merge_sccs.h"

#include "factored_transition_system.h"
#include "merge_dfp.h"

#include "../causal_graph.h"
#include "../option_parser.h"
#include "../plugin.h"
#include "../scc.h"
#include "../variable_order_finder.h"

#include <algorithm>
#include <cassert>
#include <iostream>

using namespace std;

bool compare_sccs_increasing(const unordered_set<int> &lhs, const unordered_set<int> &rhs) {
    return lhs.size() < rhs.size();
}

bool compare_sccs_decreasing(const unordered_set<int> &lhs, const unordered_set<int> &rhs) {
    return lhs.size() > rhs.size();
}

MergeSCCs::MergeSCCs(const Options &options)
    : MergeStrategy(),
      scc_order(SCCOrder(options.get_enum("scc_order"))),
      merge_order(MergeOrder(options.get_enum("merge_order"))),
      var_order_type(VariableOrderType(options.get_enum("variable_order"))),
      merge_dfp(0),
      number_of_merges_for_scc(0) {
}

MergeSCCs::~MergeSCCs() {
    delete merge_dfp;
}

void MergeSCCs::initialize(const std::shared_ptr<AbstractTask> task) {
    MergeStrategy::initialize(task);

    Options opts;
    opts.set<int>("order", 0);
    merge_dfp = new MergeDFP(opts);
    merge_dfp->initialize(task);

    TaskProxy task_proxy(*task);
    VariablesProxy vars = task_proxy.get_variables();
    int num_vars = vars.size();
    // First step: compute SCCs and put them in the right order.
    vector<vector<int>> cg;
    cg.reserve(num_vars);
    for (VariableProxy var : vars) {
        const std::vector<int> &successors =
            task_proxy.get_causal_graph().get_successors(var.get_id());
        cg.push_back(successors);
    }
    SCC scc(cg);
    const vector<vector<int>> &sccs = scc.get_result();
    if (sccs.size() == 1) {
        cout << "found single scc, continue as regular merge with the "
            "chosen external merge setting" << endl;
    } else {
        cout << "found cg sccs:" << endl;
        for (size_t i = 0; i < sccs.size(); ++i) {
            const vector<int> &single_scc = sccs[i];
            if (single_scc.size() == 1) {
                cout << "skipping scc of size 1" << endl;
            } else {
                cout << single_scc << endl;
                cg_sccs.push_back(unordered_set<int>(single_scc.begin(), single_scc.end()));
            }
        }
        switch (scc_order) {
        case TOPOLOGICAL:
            // sccs are computed in topological order
            break;
        case REVERSE_TOPOLOGICAL:
            reverse(cg_sccs.begin(), cg_sccs.end());
            break;
        case DECREASING:
            /*
              We merge starting with the *last* scc, hence sorting
              according to increasing size gives the desires decreasing
              order.
            */
            sort(cg_sccs.begin(), cg_sccs.end(), compare_sccs_increasing);
            break;
        case INCREASING:
            // see DECREASING
            sort(cg_sccs.begin(), cg_sccs.end(), compare_sccs_decreasing);
            break;
        }

        if (merge_order == DFP) {
            current_scc_ts_indices.reserve(g_variable_domain.size());
        }
    }

    // Second step: if the merge order is linear, compute the order.
    if (merge_order == LINEAR) {
        // Compute the variable order from VariableOrderFinder
        VariableOrderFinder order(task, var_order_type);
        vector<int> variable_order;
        variable_order.reserve(g_variable_domain.size());
        while (!order.done()) {
            variable_order.push_back(order.next());
        }
//        cout << "variable order finder: " << variable_order << endl;

        linear_order.reserve(g_variable_domain.size() - 1);
        int next_ts_index = g_variable_domain.size() - 1;
        unordered_map<int, int> var_to_ts_index;
        // Compute the merge order within sccs.
        while (!cg_sccs.empty()) {
            unordered_set<int> &scc = cg_sccs.back();
            int first = -1;
            int second = -1;
            int ts_index_after_scc = next_ts_index + scc.size() - 1;
            set<int> used_vars;
            for (size_t j = 0; j < variable_order.size(); ++j) {
                int var = variable_order[j];
                if (scc.count(var)) {
                    if (first == -1) {
                        first = var;
                    } else if (second == -1) {
                        second = var;
                        linear_order.push_back(make_pair(first, second));
                        ++next_ts_index;
                    } else {
                        linear_order.push_back(make_pair(next_ts_index, var));
                        ++next_ts_index;
                    }

                    // The following asserts that the internal linear order
                    // makes variables so that there is always a causal
                    // connection between the already merged variables and the
                    // next one.
                    if (!used_vars.empty()) {
                        bool connected_var = false;
                        const vector<int> &successors = cg[var];
                        for (int successor : successors) {
                            if (used_vars.count(successor)) {
                                connected_var = true;
                                break;
                            }
                        }
                        if (!connected_var) {
                            cerr << "Variable not causally connected" << endl;
                            exit_with(EXIT_CRITICAL_ERROR);
                        }
                    }
                    used_vars.insert(var);
                    scc.erase(var);
                    var_to_ts_index[var] = ts_index_after_scc;
                    if (scc.empty()) {
                        assert(ts_index_after_scc == next_ts_index);
                        break;
                    }
                }
            }
            cg_sccs.erase(cg_sccs.end());
        }
//        cout << "precomputed internal scc merge order: " << endl;
//        for (size_t i = 0; i < linear_order.size(); ++i) {
//            cout << linear_order[i].first << ", "
//                 << linear_order[i].second << endl;
//        }
//        int debug_size = linear_order.size();

        // Compute linear merge order after having merged the sccs.
        unordered_set<int> used_indices;
        int first = variable_order[0];
        if (var_to_ts_index.count(first)) {
            first = var_to_ts_index[first];
        }
        used_indices.insert(first);
        bool use_first = true;
        for (size_t i = 1; i < variable_order.size(); ++i) {
            int var = variable_order[i];
            if (var_to_ts_index.count(var)) {
                var = var_to_ts_index[var];
            }
            if (!used_indices.count(var)) {
                used_indices.insert(var);
                if (use_first) {
                    linear_order.push_back(make_pair(first, var));
                    use_first = false;
                } else {
                    linear_order.push_back(make_pair(next_ts_index, var));
                }
                ++next_ts_index;
            }
        }
//        cout << "precomputed composite scc merge order: " << endl;
//        for (size_t i = debug_size; i < linear_order.size(); ++i) {
//            cout << linear_order[i].first << ", "
//                 << linear_order[i].second << endl;
//        }
        assert(linear_order.size() == g_variable_domain.size() - 1);
    }
}

pair<int, int> MergeSCCs::get_next_dfp(
    shared_ptr<FactoredTransitionSystem> fts) {
    unordered_set<int> &current_scc = cg_sccs.back();
    pair<int, int> next_pair = merge_dfp->get_next(fts, current_scc_ts_indices);
    /*
      Try to remove both indices from the current scc. If we merge one or two
      composite transition systems resulting from previous merges of this scc,
      then no index is actually removed.
    */
    current_scc.erase(next_pair.first);
    current_scc.erase(next_pair.second);
    for (vector<int>::iterator it = current_scc_ts_indices.begin();
         it != current_scc_ts_indices.end(); ) {
      if (*it == next_pair.first || *it == next_pair.second) {
        it = current_scc_ts_indices.erase(it);
      } else {
        ++it;
      }
    }
    --remaining_merges;
    --number_of_merges_for_scc;
    return next_pair;
}

pair<int, int> MergeSCCs::get_next(
    shared_ptr<FactoredTransitionSystem> fts) {
    assert(!done());

    if (merge_order == LINEAR) {
        // use next precomputed pair
        pair<int, int> next_pair = linear_order.front();
        linear_order.erase(linear_order.begin());
        --remaining_merges;
        return next_pair;
    } else if (merge_order == DFP) {
        if (!number_of_merges_for_scc && !cg_sccs.empty()) {
            unordered_set<int> &current_scc = cg_sccs.back();
            assert(current_scc.size() > 1);
            assert(current_scc_ts_indices.empty());

            number_of_merges_for_scc = current_scc.size() - 1;
            if (number_of_merges_for_scc == 1) {
                --remaining_merges;
                --number_of_merges_for_scc;
                assert(!number_of_merges_for_scc);
                int first = *current_scc.begin();
                int second = *(++current_scc.begin());
                cg_sccs.erase(cg_sccs.end());
                return make_pair(first, second);
            }

            // Initialize current transition systems with all those contained in the scc
            int number_ts = fts->get_size();
            for (int i = 0; i < number_ts; ++i) {
                if (current_scc.count(i)) {
                    assert(fts->is_active(i));
                    current_scc_ts_indices.push_back(i);
                }
            }
            return get_next_dfp(fts);
        }

        if (number_of_merges_for_scc > 1) {
            // Add the newest transition system to the set of current ones of the scc
            current_scc_ts_indices.push_back(fts->get_size() - 1);
            return get_next_dfp(fts);
        }

        if (number_of_merges_for_scc == 1) {
            current_scc_ts_indices.push_back(fts->get_size() - 1);
            assert(current_scc_ts_indices.size() == 2);
            pair<int, int> next_pair = make_pair(current_scc_ts_indices[0],
                current_scc_ts_indices[1]);

            // Assert that we merged all transition systems that we expected to merge.
            unordered_set<int> &current_scc = cg_sccs.back();
            current_scc.erase(next_pair.first);
            current_scc.erase(next_pair.second);
            assert(current_scc.empty());
            current_scc_ts_indices.clear();
            cg_sccs.erase(cg_sccs.end());

            --remaining_merges;
            --number_of_merges_for_scc;
            assert(!number_of_merges_for_scc);
            return next_pair;
        }

        --remaining_merges;
        return merge_dfp->get_next(fts);
    } else {
        ABORT("Unknown merge order");
        return make_pair(-1, -1);
    }
}

string MergeSCCs::name() const {
    return "sccs";
}

static shared_ptr<MergeStrategy>_parse(OptionParser &parser) {
    vector<string> orders;
    orders.push_back("topological");
    orders.push_back("reverse_topological");
    orders.push_back("decreasing");
    orders.push_back("increasing");
    parser.add_enum_option("scc_order",
                           orders,
                           "choose an ordering of the sccs",
                           "topological");
    vector<string> merge_order;
    merge_order.push_back("linear");
    merge_order.push_back("dfp");
    parser.add_enum_option("merge_order",
                           merge_order,
                           "choose a merge order: linear (specify "
                           "variable_order)  or dfp.",
                           "dfp");
    // linear merge strategy option
    vector<string> variable_order;
    variable_order.push_back("CG_GOAL_LEVEL");
    variable_order.push_back("CG_GOAL_RANDOM");
    variable_order.push_back("GOAL_CG_LEVEL");
    variable_order.push_back("RANDOM");
    variable_order.push_back("LEVEL");
    variable_order.push_back("REVERSE_LEVEL");
    parser.add_enum_option("variable_order",
                           variable_order,
                           "option useful if merge_order = linear. "
                           "see VariableOrderFinder",
                           "reverse_level");
    // dfp merge strategy option
    vector<string> order;
    order.push_back("DFP");
    order.push_back("REGULAR");
    order.push_back("INVERSE");
    parser.add_enum_option("order", order, "order of transition systems", "DFP");
    Options options = parser.parse();
    if (parser.dry_run())
        return 0;
    else
        return make_shared<MergeSCCs>(options);
}

static PluginShared<MergeStrategy> _plugin("merge_sccs", _parse);
