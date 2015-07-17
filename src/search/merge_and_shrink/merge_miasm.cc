#include "merge_miasm.h"

//#include "labels.h"
//#include "transition_system.h"
//#include "merge_and_shrink_heuristic.h"
#include "sink_set_search.h"
#include "miasm_mas.h"
#include "merge_tree.h"

#include "../option_parser.h"
#include "../plugin.h"

#include <iostream>
#include <algorithm>
#include <cassert>

using namespace std;

static MergeStrategy *_parse(OptionParser &parser) {
    parser.add_option<SinkSetSearch *>(
        SinkSetSearch::option_key(),
        "search priority subsets",
        SinkSetSearch::plugin_key());

    parser.add_enum_option(MiasmInternal::option_key(),
                           MiasmInternal::S(),
                           "",
                           MiasmInternal::default_value());

    parser.add_enum_option(MiasmExternal::option_key(),
                           MiasmExternal::S(),
                           "",
                           MiasmExternal::default_value());

    Options opts = parser.parse();

    if (parser.dry_run())
        return 0;

    return new MergeMiasm(opts);
}

static Plugin<MergeStrategy> _plugin("merge_miasm", _parse);

MergeMiasm::MergeMiasm(const Options &opts)
    : MergeStrategy(),
      sink_set_search(opts.get<SinkSetSearch *>(SinkSetSearch::option_key())),
      miasm_internal(opts.get_enum(MiasmInternal::option_key())),
      miasm_external(opts.get_enum(MiasmExternal::option_key())) {
    merge_count = 0;
}

MergeMiasm::~MergeMiasm() {
}


pair<int, int> MergeMiasm::get_next(
    const vector<TransitionSystem *> &all_transition_systems) {
    /* TODO: in fact, for MIASM all work has been done
     * in the preprocess phase, including a complete merge order.
     * Thus, the transition system set is not even needed here
     * to compute a merge order. The following lines
     * (before "remaining_merge--")
     * are just for avoiding "unused variable" compiler error */
    int bound = (int)all_transition_systems.size();

    if (miasm_next[merge_count].first >= bound ||
        miasm_next[merge_count].second >= bound) {
        cerr << "impossible!";
    }

    remaining_merges--;
    return miasm_next[merge_count++];
}

void MergeMiasm::initialize(const shared_ptr<AbstractTask> task) {
    /* search for sink sets */
    sink_set_search->search();
    sink_set_search->get_sink_set(sink_sets);

    sink_set_search->miasm_abstraction->release_cache();

//    sink_set_search->reset();

//    sink_set_search->search();
//    sink_set_search->get_sink_set(sink_sets);

//    cerr << "2 before release_cache, current rss: " << getCurrentRSS() << endl;
//    sleep(5);
//    sink_set_search->miasm_abstraction->release_cache();
//    cerr << "2 the end, current rss: " << getCurrentRSS() << endl;
//    sleep(5);


    /* find the maximal weighted set packing of the priority sets */
    greedy_max_set_packing();
    cerr << "max packing" << max_packing << endl;
    /* construct the merge tree based on the max packing
     * using the internal and external merging strategy
     * specified in the options for the current MIASM */
    MiasmMergeTree miasm_tree(
        max_packing, miasm_internal, miasm_external,
        sink_set_search->get_vsir(),
        task);
    /* convert the merge tree into the merge order
     * for the convenience of of the generic MergeStrategy::get_next
     * function */
    miasm_tree.get_order(miasm_next);



//    vector<size_t> check_bounds;
//    mutex_control->get_check_bounds(miasm_next, check_bounds);
//    cerr << check_bounds << endl;
}

string MergeMiasm::name() const {
    return "miasm";
}
void MergeMiasm::dump_strategy_specific_options() const {
}

void MergeMiasm::greedy_max_set_packing() {
    max_packing.clear();
    /* the variables that have been included in the packing solution */
    set<int> included;
    /* the subsets have been sorted in the decreasing order in weight,
     * i.e., the increasing order in the R&R ratio */
    for (size_t i = 0; i < sink_sets.size(); ++i) {
        set<int> intersection;
        set_intersection(
            sink_sets[i].begin(), sink_sets[i].end(),
            included.begin(), included.end(),
            inserter(intersection, intersection.begin()));
        /* if the subset does not intersect with the set of all
         * included variables, then add this subset into the current
         * solution and add its variables as included */
        if (intersection.empty()) {
            max_packing.push_back(sink_sets[i]);
            for (set<int>::iterator j = sink_sets[i].begin();
                 j != sink_sets[i].end(); ++j) {
                included.insert(*j);
            }
        }
    }
}
