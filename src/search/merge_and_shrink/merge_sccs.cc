#include "merge_sccs.h"

#include "../causal_graph.h"
#include "../globals.h"
#include "../option_parser.h"
#include "../plugin.h"
#include "../scc.h"
#include "../variable_order_finder.h"

#include <algorithm>
#include <cassert>
#include <iostream>

using namespace std;

bool compare_sccs_increasing(const set<int> &lhs, const set<int> &rhs) {
    return lhs.size() < rhs.size();
}

bool compare_sccs_decreasing(const set<int> &lhs, const set<int> &rhs) {
    return lhs.size() > rhs.size();
}

MergeSCCs::MergeSCCs(const Options &options)
    : MergeDFP(),
      merge_order(MergeOrder(options.get_enum("merge_order"))),
      number_of_merges_for_scc(0) {

    // First step: compute SCCs and put them in the right order.
    vector<vector<int> > cg;
    cg.reserve(g_variable_domain.size());
    for (size_t var = 0; var < g_variable_domain.size(); ++var) {
        const std::vector<int> &successors = g_causal_graph->get_successors(var);
        cg.push_back(successors);
    }
    SCC scc(cg);
    const vector<vector<int> > &sccs = scc.get_result();
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
                cg_sccs.push_back(set<int>(single_scc.begin(), single_scc.end()));
            }
        }
        switch (SCCOrder(options.get_enum("scc_order"))) {
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
            current_transition_systems.reserve(g_variable_domain.size());
        }
    }

    // Second step: if the merge order is linear, compute the order.
    if (merge_order == LINEAR) {
        VariableOrderFinder order(VariableOrderType(options.get_enum("variable_order")));
        vector<int> variable_order;
        variable_order.reserve(g_variable_domain.size());
        while (!order.done()) {
            variable_order.push_back(order.next());
        }
        linear_order.reserve(g_variable_domain.size() - 1);
        int next_ts_index = g_variable_domain.size();
        unordered_map<int, int> var_to_ts_index;
        // compute the merge order within sccs.
        while (!cg_sccs.empty()) {
            set<int> &scc = cg_sccs.back();
            int first = -1;
            int second = -1;
            int ts_index_after_scc = next_ts_index + scc.size() - 1;
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
        cout << "precomputed internal scc merge order: ";
        for (size_t i = 0; i < linear_order.size(); ++i) {
            cout << linear_order[i].first << ", "
                 << linear_order[i].second << endl;
        }

        // compute linear merge order after having merged the sccs.
        int first = variable_order[0];
        if (var_to_ts_index.count(first)) {
            first = var_to_ts_index[first];
        }
        int second = variable_order[1];
        if (var_to_ts_index.count(second)) {
            second = var_to_ts_index[second];
        }
        linear_order.push_back(make_pair(first, second));
        ++next_ts_index;
        for (size_t i = 2; i < variable_order.size(); ++i) {
            int var = variable_order[i];
            if (var_to_ts_index.count(var)) {
                var = var_to_ts_index[var];
            }
            linear_order.push_back(make_pair(next_ts_index, var));
            ++next_ts_index;
        }
        assert(linear_order.size() == g_variable_domain.size() - 1);
    }
}

pair<int, int> MergeSCCs::get_next_dfp() {
    set<int> &current_scc = cg_sccs.back();
    pair<int, int> next_pair = MergeDFP::get_next(current_transition_systems);
    /*
      Try to remove both indices from the current scc. If we merge one or two
      composite transition systems resulting from previous merges of this scc,
      then no index is actually removed.
    */
    current_scc.erase(next_pair.first);
    current_scc.erase(next_pair.second);
    current_transition_systems[next_pair.first] = 0;
    current_transition_systems[next_pair.second] = 0;
    --number_of_merges_for_scc;
    return next_pair;
}

pair<int, int> MergeSCCs::get_next(const std::vector<TransitionSystem *> &all_transition_systems) {
    assert(!done());

    if (merge_order == LINEAR) {
        // use next precomputed pair
        pair<int, int> next_pair = linear_order.front();
        linear_order.erase(linear_order.begin());
        --remaining_merges;
        return next_pair;
    } else if (merge_order == DFP) {
        if (!number_of_merges_for_scc && !cg_sccs.empty()) {
            set<int> &current_scc = cg_sccs.back();
            assert(current_scc.size() > 1);
            assert(current_transition_systems.empty());

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
            for (size_t i = 0; i < all_transition_systems.size(); ++i) {
                if (current_scc.count(i)) {
                    TransitionSystem *ts = all_transition_systems[i];
                    assert(ts);
                    current_transition_systems.push_back(ts);
                } else {
                    current_transition_systems.push_back(0);
                }
            }
            return get_next_dfp();
        }

        if (number_of_merges_for_scc > 1) {
            // Add the newest transition system to the set of current ones of the scc
            current_transition_systems.push_back(all_transition_systems.back());
            return get_next_dfp();
        }

        if (number_of_merges_for_scc == 1) {
            current_transition_systems.push_back(all_transition_systems.back());
            pair<int, int> next_pair;
            bool looking_for_first = true;
            for (size_t i = 0; i < current_transition_systems.size(); ++i) {
                if (current_transition_systems[i]) {
                    if (looking_for_first) {
                        next_pair.first = i;
                        looking_for_first = false;
                    } else {
                        next_pair.second = i;
                        break;
                    }
                }
            }
            cout << "Next pair of indices: (" << next_pair.first << ", " << next_pair.second << ")" << endl;

            // Assert that we merged all transition systems that we expected to merge.
            set<int> &current_scc = cg_sccs.back();
            current_scc.erase(next_pair.first);
            current_scc.erase(next_pair.second);
            assert(current_scc.empty());
            current_transition_systems[next_pair.first] = 0;
            current_transition_systems[next_pair.second] = 0;
            for (size_t i = 0; i < current_transition_systems.size(); ++i) {
                assert(!current_transition_systems[i]);
            }
            current_transition_systems.clear();
            cg_sccs.erase(cg_sccs.end());

            --remaining_merges;
            --number_of_merges_for_scc;
            assert(!number_of_merges_for_scc);
            return next_pair;
        }

        return MergeDFP::get_next(all_transition_systems);
    } else {
        ABORT("Unknown merge order");
        return make_pair(-1, -1);
    }
}

string MergeSCCs::name() const {
    return "sccs";
}

static MergeStrategy *_parse(OptionParser &parser) {
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
    Options options = parser.parse();
    if (parser.dry_run())
        return 0;
    else
        return new MergeSCCs(options);
}

static Plugin<MergeStrategy> _plugin("merge_sccs", _parse);
