#include "merge_strategy_factory_dynamic_weighted.h"

#include "merge_dynamic_weighted.h"
#include "merge_strategy_factory_dfp.h"

#include "../task_proxy.h"

#include "../options/option_parser.h"
#include "../options/options.h"
#include "../options/plugin.h"

#include "../utils/system.h"
#include "../utils/rng.h"
#include "../utils/rng_options.h"

using namespace std;
using utils::ExitCode;

namespace merge_and_shrink {
MergeStrategyFactoryWeighted::MergeStrategyFactoryWeighted(
    const options::Options &opts)
    : MergeStrategyFactory(),
      options(opts) {
    features = utils::make_unique_ptr<Features>(opts);
}

void MergeStrategyFactoryWeighted::dump_strategy_specific_options() const {
    features->dump_weights();
}

// TODO: copied from DFP
vector<int> MergeStrategyFactoryWeighted::compute_ts_order(
    shared_ptr<AbstractTask> task,
    const options::Options &options) {
    AtomicTSOrder atomic_ts_order = AtomicTSOrder(options.get_enum("atomic_ts_order"));
    ProductTSOrder product_ts_order = ProductTSOrder(options.get_enum("product_ts_order"));
    bool atomic_before_product = options.get<bool>("atomic_before_product");
    bool randomized_order = options.get<bool>("randomized_order");
    shared_ptr<utils::RandomNumberGenerator> rng = utils::parse_rng_from_options(options);
    TaskProxy task_proxy(*task);
    int num_variables = task_proxy.get_variables().size();
    int max_transition_system_count = num_variables * 2 - 1;
    vector<int> transition_system_order;
    transition_system_order.reserve(max_transition_system_count);
    if (randomized_order) {
        for (int i = 0; i < max_transition_system_count; ++i) {
            transition_system_order.push_back(i);
        }
        rng->shuffle(transition_system_order);
    } else {
        // Compute the order in which atomic transition systems are considered
        vector<int> atomic_tso;
        for (int i = 0; i < num_variables; ++i) {
            atomic_tso.push_back(i);
        }
        if (atomic_ts_order == AtomicTSOrder::INVERSE) {
            reverse(atomic_tso.begin(), atomic_tso.end());
        } else if (atomic_ts_order == AtomicTSOrder::RANDOM) {
            rng->shuffle(atomic_tso);
        }

        // Compute the order in which product transition systems are considered
        vector<int> product_tso;
        for (int i = num_variables; i < max_transition_system_count; ++i) {
            product_tso.push_back(i);
        }
        if (product_ts_order == ProductTSOrder::NEW_TO_OLD) {
            reverse(product_tso.begin(), product_tso.end());
        } else if (product_ts_order == ProductTSOrder::RANDOM) {
            rng->shuffle(product_tso);
        }

        // Put the orders in the correct order
        if (atomic_before_product) {
            transition_system_order.insert(transition_system_order.end(),
                                           atomic_tso.begin(),
                                           atomic_tso.end());
            transition_system_order.insert(transition_system_order.end(),
                                           product_tso.begin(),
                                           product_tso.end());
        } else {
            transition_system_order.insert(transition_system_order.end(),
                                           product_tso.begin(),
                                           product_tso.end());
            transition_system_order.insert(transition_system_order.end(),
                                           atomic_tso.begin(),
                                           atomic_tso.end());
        }
    }

    return transition_system_order;
}

unique_ptr<MergeStrategy> MergeStrategyFactoryWeighted::compute_merge_strategy(
    shared_ptr<AbstractTask> task,
    FactoredTransitionSystem &fts) {
    TaskProxy task_proxy(*task);
    features->initialize(task_proxy);
    vector<int> transition_system_order = compute_ts_order(task, options);
    int max_states = options.get<int>("max_states");
    bool use_lr = options.get<bool>("use_lr");
    return utils::make_unique_ptr<MergeDynamicWeighted>(
        fts, move(features), move(transition_system_order), max_states, use_lr);
}

string MergeStrategyFactoryWeighted::name() const {
    return "dynamic merging";
}

static shared_ptr<MergeStrategyFactory>_parse(options::OptionParser &parser) {
    parser.add_option<bool>(
        "debug", "debug", "false");
    parser.add_option<int>("max_states", "shrink strategy option", "50000");
    parser.add_option<bool>("use_lr", "use label reduction", "false");

    // TS order options
    MergeStrategyFactoryDFP::add_options_to_parser(parser, false);

    // Feature weight options
    parser.add_option<int>(
        "w_causally_connected_vars",
        "prefer merging variables that are causally connected ",
        "0",
        options::Bounds("0", "100"));
    parser.add_option<int>(
        "w_bool_causally_connected_vars",
        "prefer merging variables that are causally connected ",
        "0",
        options::Bounds("0", "100"));
    parser.add_option<int>(
        "w_nonadditive_vars",
        "avoid merging additive variables",
        "0",
        options::Bounds("0", "100"));
    parser.add_option<int>(
        "w_small_transitions_states_quotient",
        "prefer merging 'sparse' transition systems",
        "0",
        options::Bounds("0", "100"));
    parser.add_option<int>(
        "w_high_transitions_states_quotient",
        "prefer merging 'dense' transition systems",
        "0",
        options::Bounds("0", "100"));
    parser.add_option<int>(
        "w_high_initial_h_value_improvement",
        "prefer merging transition systems with high initial h value improvement",
        "0",
        options::Bounds("0", "100"));
    parser.add_option<int>(
        "w_high_absolute_initial_h_value",
        "prefer merging transition systems with high absolute initial h value",
        "0",
        options::Bounds("0", "100"));
    parser.add_option<int>(
        "w_high_absolute_max_f_value",
        "prefer merging transition systems with high absolute max f value",
        "0",
        options::Bounds("0", "100"));
    parser.add_option<int>(
        "w_high_absolute_max_g_value",
        "prefer merging transition systems with high absolute max g value",
        "0",
        options::Bounds("0", "100"));
    parser.add_option<int>(
        "w_high_absolute_max_h_value",
        "prefer merging transition systems with high absolute max h value",
        "0",
        options::Bounds("0", "100"));
    parser.add_option<int>(
        "w_high_average_h_value_improvement",
        "prefer merging transition systems with high average h value",
        "0",
        options::Bounds("0", "100"));
    parser.add_option<int>(
        "w_high_initial_h_value_sum",
        "prefer merging transition systems with large number of states",
        "0",
        options::Bounds("0", "100"));
    parser.add_option<int>(
        "w_high_average_h_value_sum",
        "prefer merging transition systems with large number of edges",
        "0",
        options::Bounds("0", "100"));
    parser.add_option<int>(
        "w_dfp",
        "merge according to DFP merge strategy",
        "0",
        options::Bounds("0", "100"));
    parser.add_option<int>(
        "w_goal_relevance",
        "prefer goal relevant transition systems",
        "0",
        options::Bounds("0", "100"));
    parser.add_option<int>(
        "w_num_variables",
        "prefer transition systems with many incorporated variables",
        "0",
        options::Bounds("0", "100"));
    parser.add_option<int>(
        "w_shrink_perfectly",
        "prefer merges which allow shrinking perfectly",
        "0",
        options::Bounds("0", "100"));
    parser.add_option<int>(
        "w_num_trans",
        "prefer transition systems with few transitions",
        "0",
        options::Bounds("0", "100"));
    parser.add_option<int>(
        "w_lr_opp",
        "prefer transition systems that allow for most label reductions",
        "0",
        options::Bounds("0", "100"));
    parser.add_option<int>(
        "w_more_lr_opp",
        "prefer transition systems that allow for most label reductions",
        "0",
        options::Bounds("0", "100"));
    parser.add_option<int>(
        "w_miasm",
        "prefer transition systems that allow for most unreachable and irrelevant pruning",
        "0",
        options::Bounds("0", "100"));
    parser.add_option<int>(
        "w_mutex",
        "prefer transition systems that have facts mutex to each other",
        "0",
        options::Bounds("0", "100"));

    options::Options opts = parser.parse();
    if (opts.get<int>("w_causally_connected_vars") == 0 &&
        opts.get<int>("w_bool_causally_connected_vars") == 0 &&
        opts.get<int>("w_nonadditive_vars") == 0 &&
        opts.get<int>("w_small_transitions_states_quotient") == 0 &&
        opts.get<int>("w_high_transitions_states_quotient") == 0 &&
        opts.get<int>("w_high_initial_h_value_improvement") == 0 &&
        opts.get<int>("w_high_absolute_initial_h_value") == 0 &&
        opts.get<int>("w_high_absolute_max_f_value") == 0 &&
        opts.get<int>("w_high_absolute_max_g_value") == 0 &&
        opts.get<int>("w_high_absolute_max_h_value") == 0 &&
        opts.get<int>("w_high_average_h_value_improvement") == 0 &&
        opts.get<int>("w_high_initial_h_value_sum") == 0 &&
        opts.get<int>("w_high_average_h_value_sum") == 0 &&
        opts.get<int>("w_dfp") == 0 &&
        opts.get<int>("w_goal_relevance") == 0 &&
        opts.get<int>("w_num_variables") == 0 &&
        opts.get<int>("w_shrink_perfectly") == 0 &&
        opts.get<int>("w_num_trans") == 0 &&
        opts.get<int>("w_lr_opp") == 0 &&
        opts.get<int>("w_more_lr_opp") == 0 &&
        opts.get<int>("w_miasm") == 0 &&
        opts.get<int>("w_mutex") == 0) {
        cerr << "you must specify at least one non-zero weight!" << endl;
        utils::exit_with(ExitCode::INPUT_ERROR);
    }

    if (parser.dry_run())
        return nullptr;
    else
        return make_shared<MergeStrategyFactoryWeighted>(opts);
}

static options::PluginShared<MergeStrategyFactory> _plugin("merge_dynamic_weighted", _parse);
}
