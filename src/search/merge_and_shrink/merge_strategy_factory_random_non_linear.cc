#include "merge_strategy_factory_random_non_linear.h"

#include "merge_random_non_linear.h"

#include "../options/option_parser.h"
#include "../options/options.h"
#include "../options/plugin.h"

#include "../utils/memory.h"
#include "../utils/rng.h"
#include "../utils/rng_options.h"

#include <iostream>

using namespace std;

namespace merge_and_shrink {
MergeStrategyFactoryRandomNonLinear::MergeStrategyFactoryRandomNonLinear(
    const options::Options &options)
    : MergeStrategyFactory(),
      shrink_threshold(options.get<int>("shrink_threshold")),
      random_seed(options.get<int>("random_seed")) {
    rng = utils::parse_rng_from_options(options);
}

unique_ptr<MergeStrategy> MergeStrategyFactoryRandomNonLinear::compute_merge_strategy(
    const std::shared_ptr<AbstractTask> &,
    FactoredTransitionSystem &fts) {
    return utils::make_unique_ptr<MergeRandomNonLinear>(fts, move(rng), shrink_threshold);
}

void MergeStrategyFactoryRandomNonLinear::dump_strategy_specific_options() const {
    cout << "shrink threshold (for merging non-linearly): " << shrink_threshold << endl;
    cout << "random seed: " << random_seed << endl;
}

string MergeStrategyFactoryRandomNonLinear::name() const {
    return "random non linear";
}

static shared_ptr<MergeStrategyFactory>_parse(options::OptionParser &parser) {
    parser.document_synopsis(
        "Random non-linear merge strategy.",
        "This merge strategy randomly selects the two next transition systems "
        "among those that can be merged without shrinking (thus forcing a "
        "non-linear order on the atomic variables), before randomly "
        "merging all remaining transition systems once this is no longer "
        "possible.");
    utils::add_rng_options(parser);
    parser.add_option<int>("shrink_threshold", "shrink threshold", "50000");

    options::Options opts = parser.parse();
    if (parser.dry_run())
        return nullptr;
    else
        return make_shared<MergeStrategyFactoryRandomNonLinear>(opts);
}

static options::PluginShared<MergeStrategyFactory> _plugin("merge_random_non_linear", _parse);
}
