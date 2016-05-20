#include "merge_strategy_factory_random_linear.h"

#include "merge_random_linear.h"

#include "../task_proxy.h"

#include "../options/option_parser.h"
#include "../options/options.h"
#include "../options/plugin.h"

#include "../utils/memory.h"
#include "../utils/rng.h"
#include "../utils/rng_options.h"

#include <iostream>

using namespace std;

namespace merge_and_shrink {
MergeStrategyFactoryRandomLinear::MergeStrategyFactoryRandomLinear(
    const options::Options &options)
    : MergeStrategyFactory(),
      random_seed(options.get<int>("random_seed")) {
    rng = utils::parse_rng_from_options(options);
}

unique_ptr<MergeStrategy> MergeStrategyFactoryRandomLinear::compute_merge_strategy(
    const std::shared_ptr<AbstractTask> task,
    FactoredTransitionSystem &fts) {
    TaskProxy task_proxy(*task);
    int num_variables = task_proxy.get_variables().size();
    vector<int> randomized_variable_order(num_variables, -1);
    iota(randomized_variable_order.begin(), randomized_variable_order.end(), 0);
    rng->shuffle(randomized_variable_order);

    return utils::make_unique_ptr<MergeRandomLinear>(
        fts, move(randomized_variable_order));
}

void MergeStrategyFactoryRandomLinear::dump_strategy_specific_options() const {
    cout << "random seed: " << random_seed << endl;
}

string MergeStrategyFactoryRandomLinear::name() const {
    return "random linear";
}

static shared_ptr<MergeStrategyFactory>_parse(options::OptionParser &parser) {
    parser.document_synopsis(
        "Random linear merge strategy.",
        "This merge strategy randomly computes a variable order for merging.");
    utils::add_rng_options(parser);

    options::Options opts = parser.parse();
    if (parser.dry_run())
        return nullptr;
    else
        return make_shared<MergeStrategyFactoryRandomLinear>(opts);
}

static options::PluginShared<MergeStrategyFactory> _plugin("merge_random_linear", _parse);
}
