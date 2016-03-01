#include "merge_random_non_linear.h"

#include "factored_transition_system.h"
#include "transition_system.h"

#include "../option_parser.h"
#include "../plugin.h"

#include "../utils/memory.h"
#include "../utils/rng.h"

#include <cassert>
#include <iostream>

using namespace std;

namespace merge_and_shrink {
MergeRandomNonLinear::MergeRandomNonLinear(const Options &options)
    : MergeStrategy(),
      random_seed(options.get<int>("random_seed")),
      rng(utils::make_unique_ptr<utils::RandomNumberGenerator>(random_seed)),
      shrink_threshold(options.get<int>("shrink_threshold")) {
}

pair<int, int> MergeRandomNonLinear::get_next(
    FactoredTransitionSystem &fts) {
    assert(initialized());
    assert(!done());

    utils::RandomNumberGenerator &rng_ = *rng;

    vector<pair<int, int>> possible_noshrink_merges;
    int number_ts = fts.get_size();
    for (int ts_index1 = 0; ts_index1 < number_ts; ++ts_index1) {
        if (fts.is_active(ts_index1)) {
            for (int ts_index2 = ts_index1 + 1; ts_index2 < number_ts; ++ts_index2) {
                if (fts.is_active(ts_index2)) {
                    if (fts.get_ts(ts_index1).get_size() < shrink_threshold
                            / fts.get_ts(ts_index2).get_size()) {
                        possible_noshrink_merges.push_back(make_pair(ts_index1, ts_index2));
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
        vector<int> active_count_to_ts_index;
        for (int ts_index = 0; ts_index < number_ts; ++ts_index) {
            if (fts.is_active(ts_index)) {
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
    assert(next_index2 != -1);

    --remaining_merges;
    return make_pair(next_index1, next_index2);
}

void MergeRandomNonLinear::dump_strategy_specific_options() const {
    cout << "random seed: " << random_seed << endl;
}

string MergeRandomNonLinear::name() const {
    return "random non linear";
}

static shared_ptr<MergeStrategy>_parse(OptionParser &parser) {
    parser.document_synopsis(
        "Random non-linear merge strategy.",
        "This merge strategy randomly selects the two next transition systems "
        "among those that can be merged without shrinking (thus forcing a "
        "non-linear order on the atomic variables), before randomly "
        "merging all remaining transition systems once this is no longer "
        "possible.");
    parser.add_option<int>("random_seed", "random seed", "2015");
    parser.add_option<int>("shrink_threshold", "shrink threshold", "50000");

    Options opts = parser.parse();
    if (parser.dry_run())
        return nullptr;
    else
        return make_shared<MergeRandomNonLinear>(opts);
}

static PluginShared<MergeStrategy> _plugin("merge_random_non_linear", _parse);
}
