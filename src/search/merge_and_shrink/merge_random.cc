#include "merge_random.h"

#include "../rng.h"

#include "../option_parser.h"
#include "../plugin.h"

#include <cassert>
#include <iostream>

using namespace std;

MergeRandom::MergeRandom(const Options &options)
    : MergeStrategy(),
      random_seed(options.get<int>("random_seed")),
      rng(make_unique_ptr<RandomNumberGenerator>(random_seed)) {
}

pair<int, int> MergeRandom::get_next(const vector<TransitionSystem *> &all_transition_systems) {
    assert(initialized());
    assert(!done());
    int number_ts = all_transition_systems.size();
    vector<int> active_count_to_ts_index;
    for (int ts_index = 0; ts_index < number_ts; ++ts_index) {
        if (all_transition_systems[ts_index]) {
            active_count_to_ts_index.push_back(ts_index);
        }
    }

    RandomNumberGenerator &rng_ = *rng;
    int number_active_ts = active_count_to_ts_index.size();
    int active_index1 = rng_(number_active_ts);
    int active_index2 = rng_(number_active_ts);
    while (active_index2 == active_index1) {
        active_index2 = rng_(number_active_ts);
    }

    assert(in_bounds(active_index1, active_count_to_ts_index));
    assert(in_bounds(active_index2, active_count_to_ts_index));
    int next_index1 = active_count_to_ts_index[active_index1];
    int next_index2 = active_count_to_ts_index[active_index2];
    assert(all_transition_systems[next_index1]);
    assert(all_transition_systems[next_index2]);

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
    parser.add_option<int>("random_seed", "random seed", "2015");

    Options opts = parser.parse();
    if (parser.dry_run())
        return nullptr;
    else
        return make_shared<MergeRandom>(opts);
}

static PluginShared<MergeStrategy> _plugin("merge_random", _parse);
