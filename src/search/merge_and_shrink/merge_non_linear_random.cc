#include "merge_non_linear_random.h"

#include "transition_system.h"

#include "../rng.h"

#include "../option_parser.h"
#include "../plugin.h"

#include <cassert>
#include <iostream>

using namespace std;

MergeNonLinearRandom::MergeNonLinearRandom(const Options &options)
    : MergeStrategy(),
      random_seed(options.get<int>("random_seed")),
      rng(make_unique_ptr<RandomNumberGenerator>(random_seed)),
      shrink_threshold(options.get<int>("shrink_threshold")) {
}

pair<int, int> MergeNonLinearRandom::get_next(
    const vector<TransitionSystem *> &all_transition_systems) {
    assert(initialized());
    assert(!done());

    RandomNumberGenerator &rng_ = *rng;

    vector<pair<int, int>> possible_noshrink_merges;
    for (size_t i = 0; i < all_transition_systems.size(); ++i) {
        const TransitionSystem *ts1 = all_transition_systems[i];
        if (ts1) {
            for (size_t j = i + 1; j < all_transition_systems.size(); ++j) {
                const TransitionSystem *ts2 = all_transition_systems[j];
                if (ts2) {
                    if (ts1->get_size() < shrink_threshold / ts2->get_size()) {
                        possible_noshrink_merges.push_back(make_pair(i, j));
                    }
                }
            }
        }
    }

    int next_index1 = -1;
    int next_index2 = -1;
    if (!possible_noshrink_merges.empty()) {
        int index = rng_(possible_noshrink_merges.size());
        next_index1 = possible_noshrink_merges[index].first;
        next_index2 = possible_noshrink_merges[index].second;
    } else {
        int number_ts = all_transition_systems.size();
        vector<int> active_count_to_ts_index;
        for (int ts_index = 0; ts_index < number_ts; ++ts_index) {
            if (all_transition_systems[ts_index]) {
                active_count_to_ts_index.push_back(ts_index);
            }
        }

        int number_active_ts = active_count_to_ts_index.size();
        int active_index1 = rng_(number_active_ts);
        int active_index2 = rng_(number_active_ts);
        while (active_index2 == active_index1) {
            active_index2 = rng_(number_active_ts);
        }
        next_index1 = active_count_to_ts_index[active_index1];
        next_index2 = active_count_to_ts_index[active_index2];
    }
    assert(next_index1 != -1);
    assert(next_index1 != -2);
    assert(all_transition_systems[next_index1]);
    assert(all_transition_systems[next_index2]);

    --remaining_merges;
    cout << "Next pair of indices: (" << next_index1 << ", "
         << next_index2 << ")" << endl;
    return make_pair(next_index1, next_index2);
}

void MergeNonLinearRandom::dump_strategy_specific_options() const {
    cout << "random seed: " << random_seed << endl;
}

string MergeNonLinearRandom::name() const {
    return "non linear random";
}

static shared_ptr<MergeStrategy>_parse(OptionParser &parser) {
    parser.document_synopsis(
        "Random merge strategy.",
        "This merge strategy randomly selects the two next transition systems"
        "to merge.");
    parser.add_option<int>("random_seed", "random seed", "2015");
    parser.add_option<int>("shrink_threshold", "shrink threshold", "50000");

    Options opts = parser.parse();
    if (parser.dry_run())
        return nullptr;
    else
        return make_shared<MergeNonLinearRandom>(opts);
}

static PluginShared<MergeStrategy> _plugin("merge_non_linear_random", _parse);
