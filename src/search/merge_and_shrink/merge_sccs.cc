#include "merge_sccs.h"

#include "../causal_graph.h"
#include "../globals.h"
#include "../option_parser.h"
#include "../plugin.h"
#include "../scc.h"

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
      number_of_merges_for_scc(0) {
    vector<vector<int> > cg;
    cg.reserve(g_variable_domain.size());
    for (size_t var = 0; var < g_variable_domain.size(); ++var) {
        const std::vector<int> &successors = g_causal_graph->get_successors(var);
        cg.push_back(successors);
    }
    SCC scc(cg);
    const vector<vector<int> > &sccs = scc.get_result();
    if (sccs.size() == 1) {
        cout << "found single scc, continue as regular dfp" << endl;
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

        current_transition_systems.reserve(g_variable_domain.size() * 2 - 1);
    }
}

pair<int, int> MergeSCCs::get_next_current_scc() {
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

// TODO: split get next into linear/dfp, inherit from both base classes,
// and use according methods. merging linearly should also not be too
// complicated because we only need to merge atomic transition systems.
// as soon as all sccs have been merged, we can again choose to use further
// dfp or linear shrinking, i.e. we could have two separate options giving
// rise to four possibilities.

pair<int, int> MergeSCCs::get_next(const std::vector<TransitionSystem *> &all_transition_systems) {
    assert(!done());

    if (!number_of_merges_for_scc && !cg_sccs.empty()) {
        set<int> &current_scc = cg_sccs.back();
        assert(current_scc.size() > 1);
        assert(current_transition_systems.empty());

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

        number_of_merges_for_scc = current_scc.size() - 1;
        return get_next_current_scc();
    }

    if (number_of_merges_for_scc > 1) {
        // Add the newest transition system to the set of current ones of the scc
        current_transition_systems.push_back(all_transition_systems.back());
        return get_next_current_scc();
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
}

string MergeSCCs::name() const {
    return "sccs";
}

static MergeStrategy *_parse(OptionParser &parser) {
    vector<string> orders;
    orders.push_back("topological");
    orders.push_back("decreasing");
    orders.push_back("increasing");
    parser.add_enum_option("scc_order",
                           orders,
                           "choose an ordering of the sccs",
                           "topological");
    Options options = parser.parse();

    if (parser.dry_run())
        return 0;
    else
        return new MergeSCCs(options);
}

static Plugin<MergeStrategy> _plugin("merge_sccs", _parse);
