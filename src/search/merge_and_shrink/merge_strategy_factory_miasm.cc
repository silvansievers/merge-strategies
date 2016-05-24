#include "merge_strategy_factory_miasm.h"

#include "factored_transition_system.h"
#include "merge_miasm.h"

#include "miasm/sink_set_search.h"
#include "miasm/miasm_mas.h"

#include "../options/option_parser.h"
#include "../options/plugin.h"

#include "../utils/system.h"

#include <iostream>
#include <algorithm>
#include <cassert>

using namespace std;

namespace merge_and_shrink {
MergeStrategyFactoryMiasm::MergeStrategyFactoryMiasm(const options::Options &opts)
    : MergeStrategyFactory(),
      options(opts),
      miasm_internal(MiasmInternal(opts.get_enum("miasm_internal"))),
      miasm_external(MiasmExternal(opts.get_enum("miasm_external"))) {
    // TODO: We would like to store sink_set_search instead of options here,
    // but it requires a task object.
}

MiasmMergeTree *MergeStrategyFactoryMiasm::compute_merge_tree(
    const shared_ptr<AbstractTask> &task) {
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
    MiasmMergeTree *miasm_tree = new MiasmMergeTree(
        max_packing, miasm_internal, miasm_external,
        sink_set_search.get_vsir(),
        task);
    return miasm_tree;
}

unique_ptr<MergeStrategy> MergeStrategyFactoryMiasm::compute_merge_strategy(
    const shared_ptr<AbstractTask> &task,
    FactoredTransitionSystem &fts) {
    MiasmMergeTree *miasm_tree = compute_merge_tree(task);
    /* convert the merge tree into the merge order
     * for the convenience of of the generic MergeStrategy::get_next
     * function */
    vector<pair<int, int>> merge_order;
    miasm_tree->get_order(merge_order, TaskProxy(*task).get_variables().size());
    return utils::make_unique_ptr<MergeMiasm>(fts, move(merge_order));
}

string MergeStrategyFactoryMiasm::name() const {
    return "miasm";
}

void MergeStrategyFactoryMiasm::dump_strategy_specific_options() const {
    // TODO
}

void MergeStrategyFactoryMiasm::greedy_max_set_packing() {
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

void MergeStrategyFactoryMiasm::add_options_to_parser(options::OptionParser &parser) {
    //DEFINE_ENUM_OPT(MiasmInternal, "miasm_internal", LEVEL)
    vector<string> enum_strings;
    enum_strings.push_back("level");
    enum_strings.push_back("reverse_level");
    parser.add_enum_option("miasm_internal",
                           enum_strings,
                           "",
                           "level");

    //DEFINE_ENUM_OPT(MiasmExternal, "miasm_external", NUM_VAR_CGL)
    enum_strings.clear();
    enum_strings.push_back("num_var_cgl");
    enum_strings.push_back("rnr_size_cgl");
    enum_strings.push_back("cgrl");
    parser.add_enum_option("miasm_external",
                           enum_strings,
                           "",
                           "num_var_cgl");

    /*
      SinkSetSearch options
      This is a required option, even if it is optional for the option parser.
      This is for merge_symmetries, to avoid creating a MiasmAbstraction if
      miasm is not the fallback strategy.
    */
    parser.add_option<MiasmAbstraction *>(
        MiasmAbstraction::option_key(),
        "",
        options::OptionParser::NONE);

    //DEFINE_OPT(double, OptTimeLimit, "time_limit", "30.00")
    parser.add_option<double>("time_limit",
                              "",
                              "30.0");

    //DEFINE_OPT(size_t, OptMemoryLimit, "memory_limit", "1500000000")
    parser.add_option<int>("memory_limit",
                           "",
                           "1500000000");

    //DEFINE_OPT(int, OptSizeLimit, "size_limit", "50000")
    /** @brief An \link #DECLARE_INT_OPT int wrapper struct \endlink
     * that provides the limit on the size of an abstraction on a subset
     * that can be "enqueued" in #SinkSetSearch */
    parser.add_option<int>("size_limit",
                           "",
                           "50000");

    //DEFINE_OPT(int, OptCliqueLimit, "clique_limit", "infinity")
    parser.add_option<int>("clique_limit",
                           "",
                           "infinity");

    //DEFINE_ENUM_OPT(EnumPriority, "priority", GAIN)
    enum_strings.clear();
    enum_strings.push_back("fifo");
    enum_strings.push_back("ratio");
    enum_strings.push_back("gain");
    parser.add_enum_option("priority", enum_strings,
                           "the order in which the subsets "
                           "are dequeued in the priority queue",
                           "gain");


    //DEFINE_ENUM_OPT(EnumExpand, "expand", SINGLE)
    enum_strings.clear();
    enum_strings.push_back("single");
    enum_strings.push_back("none");
    parser.add_enum_option("expand", enum_strings,
                           "which new subsets should be added into the search"
                           "priority queue",
                           "single");

    //DEFINE_ENUM_OPT(EnumGain, "gain", ALL_ACCUR)
    enum_strings.clear();
    enum_strings.push_back("pool_guess");
    enum_strings.push_back("pool_accur");
    enum_strings.push_back("all_guess");
    enum_strings.push_back("all_accur");
    parser.add_enum_option("gain", enum_strings,
                           "",
                           "all_accur");

    //DEFINE_ENUM_OPT(EnumPrune, "prune", NONE)
    enum_strings.clear();
    enum_strings.push_back("none");
    enum_strings.push_back("cgwc_mutex");
    parser.add_enum_option("prune", enum_strings,
                           "",
                           "none");
}

static shared_ptr<MergeStrategyFactory>_parse(options::OptionParser &parser) {
    MergeStrategyFactoryMiasm::add_options_to_parser(parser);

    options::Options opts = parser.parse();

    if (parser.dry_run())
        return nullptr;

    return make_shared<MergeStrategyFactoryMiasm>(opts);
}

static options::PluginShared<MergeStrategyFactory> _plugin("merge_miasm", _parse);

}
