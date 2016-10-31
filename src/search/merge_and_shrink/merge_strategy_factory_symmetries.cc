#include "merge_strategy_factory_symmetries.h"

#include "factored_transition_system.h"
#include "merge_symmetries.h"
#include "merge_selector.h"
#include "merge_tree.h"
#include "merge_tree_factory.h"
#include "transition_system.h"

#include "symmetries/symmetry_group.h"

#include "../task_proxy.h"

#include "../options/option_parser.h"
#include "../options/options.h"
#include "../options/plugin.h"

#include "../utils/system.h"

#include <algorithm>
#include <iomanip>

using namespace std;

namespace merge_and_shrink {
MergeStrategyFactorySymmetries::MergeStrategyFactorySymmetries(
    const options::Options &options)
    : MergeStrategyFactory(),
      options(options),
      merge_tree_factory(nullptr),
      merge_selector(nullptr) {
    if (options.contains("merge_tree")) {
        merge_tree_factory = options.get<shared_ptr<MergeTreeFactory>>("merge_tree");
    }
    if (options.contains("merge_selector")) {
        merge_selector = options.get<shared_ptr<MergeSelector>>("merge_selector");
    }
}

void MergeStrategyFactorySymmetries::dump_strategy_specific_options() const {
    cout << "Options for merge symmetries:" << endl;
    cout << "    symmetries for shrinking: ";
    int symmetries_for_shrinking = options.get_enum("symmetries_for_shrinking");
    switch (symmetries_for_shrinking) {
        case 0: {
            cout << "none";
            break;
        }
        case 1: {
            cout << "atomic";
            break;
        }
        case 2: {
            cout << "local";
            break;
        }
    }
    cout << endl;
    cout << "    symmetries for merging: ";
    int symmetries_for_merging = options.get_enum("symmetries_for_merging");
    switch (symmetries_for_merging) {
        case 0: {
            cout << "none";
            break;
        }
        case 1: {
            cout << "smallest";
            break;
        }
        case 2: {
            cout << "largest";
            break;
        }
    }
    cout << endl;
    if (symmetries_for_merging) {
        cout << "    external merging: ";
        switch (options.get_enum("external_merging")) {
            case 0: {
                cout << "merge for atomic symmetry";
                break;
            }
            case 1: {
                cout << "merge for local symmetry";
                break;
            }
        }
        cout << endl;
        cout << "    internal merging: ";
        switch (options.get_enum("internal_merging")) {
            case 0: {
                cout << "linear";
                break;
            }
            case 1: {
                cout << "non linear";
                break;
            }
        }
        cout << endl;
    }
    cout << "    maxium number of m&s iterations with bliss: "
         << options.get<int>("max_bliss_iterations") << endl;
    cout << "    time limit for single bliss calls (0 means unlimited): "
         << options.get<int>("bliss_call_time_limit") << endl;
    cout << "    total time budget for bliss (0 means unlimited): "
         << options.get<int>("bliss_total_time_budget") << endl;
    cout << "    stop searching for symmetries once no symmetry was found: "
         << (options.get<bool>("stop_after_no_symmetries") ? "yes" : "no") << endl;
    cout << "    stabilize transition systems: "
         << (options.get<bool>("stabilize_transition_systems") ? "yes" : "no") << endl;
    cout << "    fallback merge strategy: ";
    if (merge_tree_factory) {
        merge_tree_factory->dump_options();
    }
    if (merge_selector) {
        merge_selector->dump_options();
    }
}

unique_ptr<MergeStrategy> MergeStrategyFactorySymmetries::compute_merge_strategy(
    const TaskProxy &task_proxy,
    FactoredTransitionSystem &fts) {
    int num_merges = task_proxy.get_variables().size() - 1;

    if (merge_selector) {
        merge_selector->initialize(task_proxy);
    }

    unique_ptr<MergeTree> merge_tree = nullptr;
    bool tree_is_miasm = false;
    if (merge_tree_factory) {
        merge_tree = merge_tree_factory->compute_merge_tree(task_proxy);
        if (merge_tree_factory->get_name() == "miasm") {
            tree_is_miasm = true;
        }
    }

    return utils::make_unique_ptr<MergeSymmetries>(
        fts,
        options,
        num_merges,
        move(merge_tree),
        merge_selector,
        tree_is_miasm);
}

string MergeStrategyFactorySymmetries::name() const {
    return "symmetries";
}

static shared_ptr<MergeStrategyFactory> _parse(options::OptionParser &parser) {
    // Options for symmetries computation
    parser.add_option<int>("max_bliss_iterations", "maximum ms iteration until "
                           "which bliss is allowed to run.",
                           "infinity");
    parser.add_option<int>("bliss_call_time_limit", "time in seconds one bliss "
                           "run is allowed to last at most (0 means no limit)",
                           "0");
    parser.add_option<int>("bliss_total_time_budget", "time in seconds bliss is "
                           "allowed to run overall (0 means no limit)",
                           "0");
    parser.add_option<bool>("stop_after_no_symmetries", "stop calling bliss "
                            "after unsuccesfull previous bliss call.",
                           "False");
    vector<string> symmetries_for_shrinking;
    symmetries_for_shrinking.push_back("NO_SHRINKING");
    symmetries_for_shrinking.push_back("ATOMIC");
    symmetries_for_shrinking.push_back("LOCAL");
    parser.add_enum_option("symmetries_for_shrinking",
                           symmetries_for_shrinking,
                           "choose the type of symmetries used for shrinking: "
                           "no shrinking, "
                           "only atomic symmetries, "
                           "local symmetries.",
                           "NO_SHRINKING");
    vector<string> symmetries_for_merging;
    symmetries_for_merging.push_back("NO_MERGING");
    symmetries_for_merging.push_back("SMALLEST");
    symmetries_for_merging.push_back("LARGEST");
    parser.add_enum_option("symmetries_for_merging",
                           symmetries_for_merging,
                           "choose the type of symmetries that should determine "
                           "the set of transition systems to be merged: "
                           "the smallest or the largest",
                           "SMALLEST");
    vector<string> external_merging;
    external_merging.push_back("MERGE_FOR_ATOMIC");
    external_merging.push_back("MERGE_FOR_LOCAL");
    parser.add_enum_option("external_merging",
                           external_merging,
                           "choose the set of transition systems to be merged: "
                           "merge for atomic: merge all transition systems affected "
                           "by the chosen symmetry, or "
                           "merge for local: merge only the transition systems "
                           "mapped (in cycles) to others. only merge every "
                           "cycle separately.",
                           "MERGE_FOR_ATOMIC");
    vector<string> internal_merging;
    internal_merging.push_back("LINEAR");
    internal_merging.push_back("NON_LINEAR");
    parser.add_enum_option("internal_merging",
                           internal_merging,
                           "choose the order in which to merge the set of "
                           "transition systems to be merged (only useful with "
                           "MERGE_FOR_ATOMIC): "
                           "linear (obvious), "
                           "non linear, which means to first merge every cycle, "
                           "and then the resulting intermediate transition systems.",
                           "LINEAR");

    // Options for GraphCreator
    parser.add_option<bool>("stabilize_transition_systems", "compute symmetries that "
                            "stabilize transition systems, i.e. that are local.", "false");
    parser.add_option<bool>("debug_graph_creator", "produce dot readable output "
                            "from the graph generating methods", "false");

    // Fallback strategy options
    parser.add_option<shared_ptr<MergeTreeFactory>>(
        "merge_tree",
        "the fallback merge 'strategy' to use if a precomputed strategy should"
        "be used.",
        options::OptionParser::NONE);
    parser.add_option<shared_ptr<MergeSelector>>(
        "merge_selector",
        "the fallback merge 'strategy' to use if a stateless strategy should"
        "be used.",
        options::OptionParser::NONE);

    options::Options options = parser.parse();
    if (options.get<int>("bliss_call_time_limit")
            && options.get<int>("bliss_total_time_budget")) {
        cerr << "Please only specify bliss_call_time_limit or "
                "bliss_total_time_budget but not both" << endl;
        utils::exit_with(utils::ExitCode::CRITICAL_ERROR);
    }
    if (options.get_enum("symmetries_for_shrinking") == 0
            && options.get_enum("symmetries_for_merging") == 0) {
        cerr << "Please use symmetries at least for shrinking or merging." << endl;
        utils::exit_with(utils::ExitCode::CRITICAL_ERROR);
    }
    bool merge_tree = options.contains("merge_tree");
    bool merge_selector = options.contains("merge_selector");
    if ((merge_tree && merge_selector) || (!merge_tree && !merge_selector)) {
        cerr << "You have to specify exactly one of the options merge_tree "
                "and merg_selector!" << endl;
        utils::exit_with(utils::ExitCode::INPUT_ERROR);
    }
    if (parser.dry_run())
        return nullptr;
    else
        return make_shared<MergeStrategyFactorySymmetries>(options);
}

static options::PluginShared<MergeStrategyFactory> _plugin("merge_symmetries", _parse);
}
