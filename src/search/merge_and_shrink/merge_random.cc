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
    int size = all_transition_systems.size();
    RandomNumberGenerator &rng_ = *rng;

    int index1 = -1;
    while (true) {
        index1 = rng_(size);
        if (all_transition_systems[index1]) {
            break;
        }
    }

    int index2 = -1;
    while (true) {
        index2 = rng_(size);
        if (all_transition_systems[index2] && index2 != index1) {
            break;
        }
    }

    --remaining_merges;
    return make_pair(index1, index2);
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
