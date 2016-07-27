#include "merge_strategy_factory_sccs.h"

#include "factored_transition_system.h"
#include "merge_sccs.h"
#include "merge_scoring_function_dfp.h"
#include "merge_scoring_function_goal_relevance.h"
#include "merge_scoring_function_single_random.h"
#include "merge_scoring_function_total_order.h"
#include "merge_selector_score_based_filtering.h"
#include "merge_tree_factory_linear.h"
#include "transition_system.h"

#include "../causal_graph.h"
#include "../scc.h"
#include "../variable_order_finder.h"

#include "../options/option_parser.h"
#include "../options/plugin.h"

#include "../utils/logging.h"

#include <algorithm>
#include <cassert>
#include <iostream>

using namespace std;

namespace merge_and_shrink {
bool compare_sccs_increasing(const vector<int> &lhs, const vector<int> &rhs) {
    return lhs.size() < rhs.size();
}

bool compare_sccs_decreasing(const vector<int> &lhs, const vector<int> &rhs) {
    return lhs.size() > rhs.size();
}

MergeStrategyFactorySCCs::MergeStrategyFactorySCCs(const options::Options &options)
    : MergeStrategyFactory(), options(options) {
}

unique_ptr<MergeStrategy> MergeStrategyFactorySCCs::compute_merge_strategy(
        shared_ptr<AbstractTask> task,
        FactoredTransitionSystem &fts) {
    OrderOfSCCs order_of_sccs(static_cast<OrderOfSCCs>(options.get_enum("order_of_sccs")));
    InternalMergeOrder internal_merge_order(static_cast<InternalMergeOrder>(options.get_enum("internal_merge_order")));
    vector<int> linear_variable_order;
    vector<vector<int>> non_singleton_cg_sccs;
    vector<int> indices_of_merged_sccs;

    TaskProxy task_proxy(*task);
    VariablesProxy vars = task_proxy.get_variables();
    int num_vars = vars.size();

    shared_ptr<MergeSelectorScoreBasedFiltering> dfp_selector = nullptr;
    if (internal_merge_order == InternalMergeOrder::DFP) {
        vector<shared_ptr<MergeScoringFunction>> scoring_functions;
        scoring_functions.push_back(make_shared<MergeScoringFunctionGoalRelevance>());
        scoring_functions.push_back(make_shared<MergeScoringFunctionDFP>());

        bool randomized_order = options.get<bool>("randomized_order");
        if (randomized_order) {
            shared_ptr<MergeScoringFunctionSingleRandom> scoring_random =
                make_shared<MergeScoringFunctionSingleRandom>(options);
            scoring_functions.push_back(scoring_random);
        } else {
            shared_ptr<MergeScoringFunctionTotalOrder> scoring_total_order =
                make_shared<MergeScoringFunctionTotalOrder>(options);
            scoring_functions.push_back(scoring_total_order);
        }
        dfp_selector = make_shared<MergeSelectorScoreBasedFiltering>(
            move(scoring_functions));
        dfp_selector->initialize(task);
    }
    if (internal_merge_order == InternalMergeOrder::LINEAR) {
        VariableOrderFinder vof(task, VariableOrderType(options.get_enum("variable_order")));
        linear_variable_order.reserve(num_vars);
        while (!vof.done()) {
            linear_variable_order.push_back(vof.next());
        }
        cout << "linear variable order: " << linear_variable_order << endl;
    }

    // Compute SCCs of the causal graph
    vector<vector<int>> cg;
    cg.reserve(num_vars);
    for (VariableProxy var : vars) {
        const vector<int> &successors =
            task_proxy.get_causal_graph().get_successors(var.get_id());
        cg.push_back(successors);
    }
//    cout << "CG:" << endl;
//    for (size_t var = 0; var < cg.size(); ++var) {
//        cout << var << ": " << cg[var] << endl;
//    }
    SCC scc(cg);
    vector<vector<int>> sccs(scc.get_result());

    // Put the SCCs in the desired order
    switch (order_of_sccs) {
    case OrderOfSCCs::TOPOLOGICAL:
        // SCCs are computed in topological order
        break;
    case OrderOfSCCs::REVERSE_TOPOLOGICAL:
        // SCCs are computed in topological order
        reverse(sccs.begin(), sccs.end());
        break;
    case OrderOfSCCs::DECREASING:
        sort(sccs.begin(), sccs.end(), compare_sccs_decreasing);
        break;
    case OrderOfSCCs::INCREASING:
        sort(sccs.begin(), sccs.end(), compare_sccs_increasing);
        break;
    }

    /*
      Compute the indices at which the merged SCCs can be found when all
      SCCs have been merged.
    */
    int index = num_vars - 1;
    cout << "found cg sccs:" << endl;
    indices_of_merged_sccs.reserve(sccs.size());
    for (const vector<int> &scc : sccs) {
        cout << scc << endl;
        int scc_size = scc.size();
        if (scc_size == 1) {
            indices_of_merged_sccs.push_back(scc.front());
        } else {
            index += scc_size - 1;
            indices_of_merged_sccs.push_back(index);
            // only store non-singleton sccs for internal merging
            non_singleton_cg_sccs.push_back(vector<int>(scc.begin(), scc.end()));
        }
    }
    if (sccs.size() == 1) {
        cout << "Only one single SCC" << endl;
    }
    if (static_cast<int>(sccs.size()) == num_vars) {
        cout << "Only singleton SCCs" << endl;
        assert(non_singleton_cg_sccs.empty());
    }
//    cout << "indices of merged sccs: " << indices_of_merged_sccs << endl;

    return utils::make_unique_ptr<MergeSCCs>(
        fts,
        internal_merge_order,
        move(linear_variable_order),
        dfp_selector,
        move(non_singleton_cg_sccs),
        move(indices_of_merged_sccs));
}

void MergeStrategyFactorySCCs::dump_strategy_specific_options() const {
    OrderOfSCCs order_of_sccs(static_cast<OrderOfSCCs>(options.get_enum("order_of_sccs")));
    cout << "Merge order of sccs: ";
    switch (order_of_sccs) {
    case OrderOfSCCs::TOPOLOGICAL:
        cout << "topological";
        break;
    case OrderOfSCCs::REVERSE_TOPOLOGICAL:
        cout << "reverse topological";
        break;
    case OrderOfSCCs::DECREASING:
        cout << "decreasing";
        break;
    case OrderOfSCCs::INCREASING:
        cout << "increasing";
        break;
    }
    cout << endl;

    InternalMergeOrder internal_merge_order(static_cast<InternalMergeOrder>(
        options.get_enum("internal_merge_order")));
    cout << "Internal merge order for sccs: ";
    switch (internal_merge_order) {
    case InternalMergeOrder::LINEAR:
        cout << "linear";
        break;
    case InternalMergeOrder::DFP:
        cout << "dfp";
        break;
    }
    cout << endl;
}

string MergeStrategyFactorySCCs::name() const {
    return "sccs";
}

static shared_ptr<MergeStrategyFactory>_parse(options::OptionParser &parser) {
    vector<string> order_of_sccs;
    order_of_sccs.push_back("topological");
    order_of_sccs.push_back("reverse_topological");
    order_of_sccs.push_back("decreasing");
    order_of_sccs.push_back("increasing");
    parser.add_enum_option("order_of_sccs",
                           order_of_sccs,
                           "choose an ordering of the sccs: topological, (cg "
                           "order) reverse_topological (cg order), decreasing "
                           "or increasing (in size).");
    vector<string> internal_merge_order;
    internal_merge_order.push_back("linear");
    internal_merge_order.push_back("dfp");
    parser.add_enum_option("internal_merge_order",
                           internal_merge_order,
                           "choose an internal merge order: linear (specify "
                           "variable_order) or dfp (specify dfp order options).",
                           "dfp");

    // linear
    MergeTreeFactoryLinear::add_options_to_parser(parser);

    // dfp
    MergeScoringFunctionTotalOrder::add_options_to_parser(parser);
    parser.add_option<bool>(
        "randomized_order",
        "If true, use a 'globally' randomized order, i.e. all transition "
        "systems are considered in an arbitrary order. This renders all other "
        "ordering options void.",
        "false");

    options::Options options = parser.parse();
    if (parser.dry_run())
        return 0;
    else
        return make_shared<MergeStrategyFactorySCCs>(options);
}

static options::PluginShared<MergeStrategyFactory> _plugin("merge_sccs", _parse);
}
