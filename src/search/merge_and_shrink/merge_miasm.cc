#include "merge_miasm.h"

#include "factored_transition_system.h"

#include "miasm/sink_set_search.h"
#include "miasm/miasm_mas.h"
#include "miasm/merge_tree.h"

#include "../plugin.h"

#include <iostream>
#include <algorithm>
#include <cassert>

using namespace std;

DEFINE_OPT(double, OptTimeLimit, "time_limit", "30.00")
DEFINE_OPT(size_t, OptMemoryLimit, "memory_limit", "1500000000")
DEFINE_OPT(int, OptSizeLimit, "size_limit", "50000")
DEFINE_OPT(int, OptCliqueLimit, "clique_limit", "infinity")

#define X(a) # a

DEFINE_ENUM_OPT(EnumPriority, "priority", GAIN)

DEFINE_ENUM_OPT(EnumExpand, "expand", SINGLE)

DEFINE_ENUM_OPT(EnumGain, "gain", ALL_ACCUR)

DEFINE_ENUM_OPT(EnumPrune, "prune", NONE)

#undef X

namespace merge_and_shrink {
static shared_ptr<MergeStrategy>_parse(OptionParser &parser) {
    parser.add_enum_option(MiasmInternal::option_key(),
                           MiasmInternal::S(),
                           "",
                           MiasmInternal::default_value());

    parser.add_enum_option(MiasmExternal::option_key(),
                           MiasmExternal::S(),
                           "",
                           MiasmExternal::default_value());

    //SinkSetSearch options
    parser.add_option<MiasmAbstraction *>(MiasmAbstraction::option_key(),
                                          "",
                                          MiasmAbstraction::plugin_key());

    parser.add_option<double>(OptTimeLimit::opt_key(),
                              "",
                              OptTimeLimit::def_val());

    parser.add_option<int>(OptMemoryLimit::opt_key(),
                           "",
                           OptMemoryLimit::def_val());

    parser.add_option<int>(OptSizeLimit::opt_key(),
                           "",
                           OptSizeLimit::def_val());

    parser.add_option<int>(OptCliqueLimit::opt_key(),
                           "",
                           OptCliqueLimit::def_val());

    parser.add_enum_option(EnumPriority::option_key(), EnumPriority::S(),
                           "the order in which the subsets "
                           "are dequeued in the priority queue",
                           EnumPriority::default_value());

    parser.add_enum_option(EnumExpand::option_key(), EnumExpand::S(),
                           "which new subsets should be added into the search"
                           "priority queue",
                           EnumExpand::default_value());


    parser.add_enum_option(EnumGain::option_key(), EnumGain::S(),
                           "",
                           EnumGain::default_value());

    parser.add_enum_option(EnumPrune::option_key(), EnumPrune::S(),
                           "",
                           EnumPrune::default_value());

    Options opts = parser.parse();

    if (parser.dry_run())
        return nullptr;

    return make_shared<MergeMiasm>(opts);
}

static PluginShared<MergeStrategy> _plugin("merge_miasm", _parse);

MergeMiasm::MergeMiasm(const Options &opts)
    : MergeStrategy(),
      options(opts),
      miasm_internal(opts.get_enum(MiasmInternal::option_key())),
      miasm_external(opts.get_enum(MiasmExternal::option_key())) {
    merge_count = 0;
}

MergeMiasm::~MergeMiasm() {
}


pair<int, int> MergeMiasm::get_next(FactoredTransitionSystem &fts) {
    /* TODO: in fact, for MIASM all work has been done
     * in the preprocess phase, including a complete merge order.
     * Thus, the transition system set is not even needed here
     * to compute a merge order. The following lines
     * (before "remaining_merge--")
     * are just for avoiding "unused variable" compiler error */
    if (!fts.is_active(miasm_next[merge_count].first) ||
        !fts.is_active(miasm_next[merge_count].second)) {
        ABORT("Invalid next merge index");
    }

    remaining_merges--;
    return miasm_next[merge_count++];
}

void MergeMiasm::initialize(const shared_ptr<AbstractTask> task) {
    MergeStrategy::initialize(task);
    /* search for sink sets */
    SinkSetSearch sink_set_search(options, task);
    sink_set_search.search();
    sink_set_search.get_sink_set(sink_sets);

    sink_set_search.miasm_abstraction->release_cache();

    /* find the maximal weighted set packing of the priority sets */
    greedy_max_set_packing();
//    cerr << "max packing" << max_packing << endl;
    /* construct the merge tree based on the max packing
     * using the internal and external merging strategy
     * specified in the options for the current MIASM */
    MiasmMergeTree miasm_tree(
        max_packing, miasm_internal, miasm_external,
        sink_set_search.get_vsir(),
        task);
    /* convert the merge tree into the merge order
     * for the convenience of of the generic MergeStrategy::get_next
     * function */
    miasm_tree.get_order(miasm_next, TaskProxy(*task).get_variables().size());
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
}
