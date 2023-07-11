#include "merge_tree_factory_miasm.h"

#include "merge_selector.h"
#include "merge_tree.h"

#include "miasm/sink_set_search.h"
#include "miasm/merge_tree.h"
#include "miasm/miasm_mas.h"

#include "../plugins/plugin.h"

#include "../utils/system.h"

#include <iostream>
#include <algorithm>
#include <cassert>

using namespace std;

namespace merge_and_shrink {
MergeTreeFactoryMiasm::MergeTreeFactoryMiasm(const plugins::Options &opts)
    : MergeTreeFactory(opts),
      options(opts),
      miasm_internal(opts.get<MiasmInternal>("miasm_internal")),
      miasm_external(opts.get<MiasmExternal>("miasm_external")),
      fallback_merge_selector(nullptr),
      trivial_partitioning(false),
      log(utils::get_log_from_options(opts)) {
    // TODO: We would like to store sink_set_search instead of options here,
    // but it requires a task object.
    if (options.contains("fallback_merge_selector")) {
        fallback_merge_selector = options.get<shared_ptr<MergeSelector>>("fallback_merge_selector");
    }
}

MiasmMergeTree *MergeTreeFactoryMiasm::compute_miasm_merge_tree(
    const TaskProxy &task_proxy) {
    /* search for sink sets */
    SinkSetSearch sink_set_search(options, task_proxy);
    sink_set_search.search();
    sink_set_search.get_sink_set(sink_sets);

    sink_set_search.miasm_abstraction->release_cache();

    /* find the maximal weighted set packing of the priority sets */
    greedy_max_set_packing();
//    cerr << "max packing" << max_packing << endl;
    if (max_packing.size() == task_proxy.get_variables().size()) {
        trivial_partitioning = true;
        log << "Found a trivial variable partitioning, ";
        if (fallback_merge_selector) {
            log << "using fallback merge strategy" << endl;
        } else {
            log << "but no fallback merge strategy specified." << endl;
        }
    }

    /* construct the merge tree based on the max packing
     * using the internal and external merging strategy
     * specified in the options for the current MIASM */
    MiasmMergeTree *miasm_tree = new MiasmMergeTree(
        max_packing, miasm_internal, miasm_external,
        sink_set_search.get_vsir(),
        task_proxy);
    return miasm_tree;
}

unique_ptr<MergeTree> MergeTreeFactoryMiasm::compute_merge_tree(
    const TaskProxy &task_proxy) {
    int num_ts = task_proxy.get_variables().size();

    // compute the merge tree in MiasmMergeTree form
    MiasmMergeTree *miasm_tree = compute_miasm_merge_tree(task_proxy);

    // get the actual merge order
    vector<pair<int, int>> merge_order;
    int next_ts_index = num_ts;
    while (true) {
        pair<int, int> next_merge = miasm_tree->select_next_and_update(next_ts_index);
        if (next_merge.first == -1) {
            break;
        }
        merge_order.push_back(next_merge);
        ++next_ts_index;
    }
    delete miasm_tree;

    // compute the merge tree in MergeTree form from the order
    // TODO: change the miasm computation to use it directly!
    map<int, MergeTreeNode *> index_to_tree;
    for (int atomic_ts_index = 0; atomic_ts_index < num_ts; ++atomic_ts_index) {
        index_to_tree[atomic_ts_index] = new MergeTreeNode(atomic_ts_index);
    }
    next_ts_index = num_ts;
    for (const pair<int, int> &merge : merge_order) {
        int ts_index1 = merge.first;
        int ts_index2 = merge.second;
        index_to_tree[next_ts_index] =
            new MergeTreeNode(index_to_tree[ts_index1], index_to_tree[ts_index2]);
        ++next_ts_index;
    }
    MergeTreeNode *root = index_to_tree[next_ts_index - 1];
//    MergeTree merge_tree(root, g_rng());
//    vector<pair<int, int>> other_merge_order;
//    next_ts_index = num_ts;
//    while (!merge_tree.done()) {
//        other_merge_order.push_back(merge_tree.get_next_merge(next_ts_index));
//        ++next_ts_index;
//    }
//    assert(merge_order == other_merge_order);

    return utils::make_unique_ptr<MergeTree>(root, rng, update_option);
}

string MergeTreeFactoryMiasm::name() const {
    return "miasm";
}

void MergeTreeFactoryMiasm::dump_tree_specific_options(utils::LogProxy &) const {
    // TODO
}

void MergeTreeFactoryMiasm::greedy_max_set_packing() {
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

void MergeTreeFactoryMiasm::add_options_to_feature(plugins::Feature &feature) {
    MergeTreeFactory::add_options_to_feature(feature);

    feature.add_option<MiasmInternal>("miasm_internal",
                                          "",
                                          "level");

    feature.add_option<MiasmExternal>("miasm_external",
                                          "",
                                          "num_var_cgl");

    /*
      SinkSetSearch options
      This is a required option, even if it is optional for the option parser.
      This is for merge_symmetries, to avoid creating a MiasmAbstraction if
      miasm is not the fallback strategy.
    */
    feature.add_option<shared_ptr<MiasmAbstraction>>(
        MiasmAbstraction::option_key(),
        "",
        plugins::ArgumentInfo::NO_DEFAULT);

    //DEFINE_OPT(double, OptTimeLimit, "time_limit", "30.00")
    feature.add_option<double>("time_limit",
                              "",
                              "30.0");

    //DEFINE_OPT(size_t, OptMemoryLimit, "memory_limit", "1500000000")
    feature.add_option<int>("memory_limit",
                           "",
                           "1500000000");

    //DEFINE_OPT(int, OptSizeLimit, "size_limit", "50000")
    /** @brief An \link #DECLARE_INT_OPT int wrapper struct \endlink
     * that provides the limit on the size of an abstraction on a subset
     * that can be "enqueued" in #SinkSetSearch */
    feature.add_option<int>("size_limit",
                           "",
                           "50000");

    //DEFINE_OPT(int, OptCliqueLimit, "clique_limit", "infinity")
    feature.add_option<int>("clique_limit",
                           "",
                           "infinity");

    feature.add_option<EnumPriority>("priority",
                                         "the order in which the subsets "
                                         "are dequeued in the priority queue",
                                         "gain");


    feature.add_option<EnumExpand>("expand",
                                       "which new subsets should be added into the search"
                                       "priority queue",
                                       "single");

    feature.add_option<EnumGain>("gain",
                                     "",
                                     "all_accur");

    feature.add_option<EnumPrune>("prune",
                                      "",
                                      "none");

    feature.add_option<shared_ptr<MergeSelector>>(
        "fallback_merge_selector",
        "the fallback merge 'strategy' to use if a stateless strategy should"
        "be used.",
        plugins::ArgumentInfo::NO_DEFAULT);

    utils::add_log_options_to_feature(feature);
}

class MergeTreeFactoryMiasmFeature : public plugins::TypedFeature<MergeTreeFactory, MergeTreeFactoryMiasm> {
public:
    MergeTreeFactoryMiasmFeature() : TypedFeature("miasm") {
        MergeTreeFactory::add_options_to_feature(*this);
        MergeTreeFactoryMiasm::add_options_to_feature(*this);
    }
};

static plugins::FeaturePlugin<MergeTreeFactoryMiasmFeature> _plugin;

//DEFINE_ENUM_OPT(MiasmInternal, "miasm_internal", LEVEL)
static plugins::TypedEnumPlugin<MiasmInternal> _enum_miasminternal_plugin({
      {"level", ""},
      {"reverse_level", ""},
  });

//DEFINE_ENUM_OPT(MiasmExternal, "miasm_external", NUM_VAR_CGL)
static plugins::TypedEnumPlugin<MiasmExternal> _enum_miasmexternal_plugin({
      {"num_var_cgl", ""},
      {"rnr_size_cgl", ""},
      {"cgrl", ""},
  });

//DEFINE_ENUM_OPT(EnumPriority, "priority", GAIN)
static plugins::TypedEnumPlugin<EnumPriority> _enum_priority_plugin({
      {"fifo", ""},
      {"ratio", ""},
      {"gain", ""},
  });

//DEFINE_ENUM_OPT(EnumExpand, "expand", SINGLE)
static plugins::TypedEnumPlugin<EnumExpand> _enum_expand_plugin({
      {"single", ""},
      {"none", ""},
  });

//DEFINE_ENUM_OPT(EnumGain, "gain", ALL_ACCUR)
static plugins::TypedEnumPlugin<EnumGain> _enum_gain_plugin({
    {"pool_guess", ""},
    {"pool_accur", ""},
    {"all_guess", ""},
    {"all_accur", ""},
});

//DEFINE_ENUM_OPT(EnumPrune, "prune", NONE)
static plugins::TypedEnumPlugin<EnumPrune> _enum_prune_plugin({
    {"none", ""},
    {"cgwc_mutex", ""},
});
}
