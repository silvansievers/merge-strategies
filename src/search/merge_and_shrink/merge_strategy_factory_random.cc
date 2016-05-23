#include "merge_strategy_factory_random.h"

#include "merge_random.h"

#include "../options/option_parser.h"
#include "../options/options.h"
#include "../options/plugin.h"

#include "../utils/memory.h"
#include "../utils/rng.h"
#include "../utils/rng_options.h"

#include <iostream>

using namespace std;

namespace merge_and_shrink {
MergeStrategyFactoryRandom::MergeStrategyFactoryRandom(
    const options::Options &options)
    : MergeStrategyFactory(),
      random_seed(options.get<int>("random_seed")) {
    rng = utils::parse_rng_from_options(options);
}

unique_ptr<MergeStrategy> MergeStrategyFactoryRandom::compute_merge_strategy(
    const std::shared_ptr<AbstractTask> &,
    FactoredTransitionSystem &fts) {
    return utils::make_unique_ptr<MergeRandom>(fts, move(rng));
}

void MergeStrategyFactoryRandom::dump_strategy_specific_options() const {
    cout << "random seed: " << random_seed << endl;
}

string MergeStrategyFactoryRandom::name() const {
    return "random";
}

static shared_ptr<MergeStrategyFactory>_parse(options::OptionParser &parser) {
    parser.document_synopsis(
        "Random merge strategy.",
        "This merge strategy randomly selects the two next transition systems"
        "to merge.");
    utils::add_rng_options(parser);

    options::Options opts = parser.parse();
    if (parser.dry_run())
        return nullptr;
    else
        return make_shared<MergeStrategyFactoryRandom>(opts);
}

static options::PluginShared<MergeStrategyFactory> _plugin("merge_random", _parse);
}
