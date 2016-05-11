#include "merge_random.h"

#include "factored_transition_system.h"

#include "../option_parser.h"
#include "../plugin.h"

#include "../utils/collections.h"
#include "../utils/memory.h"
#include "../utils/rng.h"
#include "../utils/rng_options.h"

#include <cassert>
#include <iostream>

using namespace std;

namespace merge_and_shrink {
MergeRandom::MergeRandom(const Options &options)
    : MergeStrategy(),
      random_seed(options.get<int>("random_seed")) {
    rng = utils::parse_rng_from_options(options);
}

pair<int, int> MergeRandom::get_next(
    FactoredTransitionSystem &fts) {
    assert(initialized());
    assert(!done());
    int number_ts = fts.get_size();
    vector<int> active_count_to_ts_index;
    for (int ts_index = 0; ts_index < number_ts; ++ts_index) {
        if (fts.is_active(ts_index)) {
            active_count_to_ts_index.push_back(ts_index);
        }
    }

    utils::RandomNumberGenerator &rng_ = *rng;
    int number_active_ts = active_count_to_ts_index.size();
    int active_index1 = rng_(number_active_ts);
    int active_index2 = rng_(number_active_ts);
    while (active_index2 == active_index1) {
        active_index2 = rng_(number_active_ts);
    }

    assert(utils::in_bounds(active_index1, active_count_to_ts_index));
    assert(utils::in_bounds(active_index2, active_count_to_ts_index));
    int next_index1 = active_count_to_ts_index[active_index1];
    int next_index2 = active_count_to_ts_index[active_index2];

    --remaining_merges;
    return make_pair(next_index1, next_index2);
}

void MergeRandom::dump_strategy_specific_options() const {
    cout << "random seed: " << random_seed << endl;
}

string MergeRandom::name() const {
    return "random";
}

static shared_ptr<MergeStrategy>_parse(OptionParser &parser) {
    parser.document_synopsis(
        "Random merge strategy.",
        "This merge strategy randomly selects the two next transition systems"
        "to merge.");
    utils::add_rng_options(parser);

    Options opts = parser.parse();
    if (parser.dry_run())
        return nullptr;
    else
        return make_shared<MergeRandom>(opts);
}

static PluginShared<MergeStrategy> _plugin("merge_random", _parse);
}
