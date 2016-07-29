#include "merge_scoring_function_miasm.h"

#include "factored_transition_system.h"
#include "shrink_bisimulation.h"
#include "transition_system.h"

#include "../options/option_parser.h"
#include "../options/options.h"
#include "../options/plugin.h"

#include "../utils/logging.h"
#include "../utils/math.h"

#include <math.h>

using namespace std;

namespace merge_and_shrink {
MergeScoringFunctionMIASM::MergeScoringFunctionMIASM(
    const options::Options &options)
    : max_states(options.get<int>("max_states")) {
}

// TODO: copied from MergeAndShrinkHeuristic
pair<int, int> MergeScoringFunctionMIASM::compute_shrink_sizes(
    int size1, int size2) const {
    // Bound both sizes by max allowed size before merge.
    int new_size1 = min(size1, max_states);
    int new_size2 = min(size2, max_states);

    if (!utils::is_product_within_limit(new_size1, new_size2, max_states)) {
        int balanced_size = int(sqrt(max_states));

        if (new_size1 <= balanced_size) {
            // Size of the first transition system is small enough. Use whatever
            // is left for the second transition system.
            new_size2 = max_states / new_size1;
        } else if (new_size2 <= balanced_size) {
            // Inverted case as before.
            new_size1 = max_states / new_size2;
        } else {
            // Both transition systems are too big. We set both target sizes
            // to balanced_size. An alternative would be to set one to
            // N1 = balanced_size and the other to N2 = max_states /
            // balanced_size, to get closer to the allowed maximum.
            // However, this would make little difference (N2 would
            // always be N1, N1 + 1 or N1 + 2), and our solution has the
            // advantage of treating the transition systems symmetrically.
            new_size1 = balanced_size;
            new_size2 = balanced_size;
        }
    }
    assert(new_size1 <= size1 && new_size2 <= size2);
    assert(new_size1 <= max_states);
    assert(new_size2 <= max_states);
    assert(new_size1 * new_size2 <= max_states);
    return make_pair(new_size1, new_size2);
}

vector<double> MergeScoringFunctionMIASM::compute_scores(
    FactoredTransitionSystem &fts,
    const vector<pair<int, int>> &merge_candidates) {
    vector<double> scores;
    scores.reserve(merge_candidates.size());
    for (pair<int, int> merge_candidate : merge_candidates ) {
        int ts_index1 = merge_candidate.first;
        int ts_index2 = merge_candidate.second;

        // Copy the transition systems (distances etc)
        int copy_ts_index1 = fts.copy(ts_index1);
        int copy_ts_index2 = fts.copy(ts_index2);
        pair<int, int> shrink_sizes =
            compute_shrink_sizes(fts.get_ts(copy_ts_index1).get_size(),
                                 fts.get_ts(copy_ts_index2).get_size());

        // shrink before merge (with implicit threshold = 1,
        // i.e. always try to shrink)
        options::Options options;
        options.set<bool>("greedy", false);
        options.set<int>("at_limit", 0);
        bool silent = true;
        ShrinkBisimulation shrink_bisim(options);
        shrink_bisim.shrink(fts, copy_ts_index1, shrink_sizes.first, silent);
        shrink_bisim.shrink(fts, copy_ts_index2, shrink_sizes.second, silent);

        // perform the merge
        int merge_index = fts.merge(copy_ts_index1, copy_ts_index2, true, false);

        // return 0 if the merge is unsolvable (i.e. empty transition system)
        double score = 0;
        if (fts.get_ts(merge_index).is_solvable()) {
            int expected_size = fts.get_ts(ts_index1).get_size() *
                fts.get_ts(ts_index2).get_size();
            assert(expected_size);
            int new_size = fts.get_ts(merge_index).get_size();
            assert(new_size <= expected_size);
            score = static_cast<double>(new_size) /
                static_cast<double>(expected_size);
        }
        scores.push_back(score);

        // delete the merge and reset
        fts.release_copies();
    }
    return scores;
}

string MergeScoringFunctionMIASM::name() const {
    return "miasm";
}

static shared_ptr<MergeScoringFunction>_parse(options::OptionParser &parser) {
    parser.document_synopsis(
        "miasm",
        "This scoring function for every candidate merge computes the "
        "actual merge and computes the ratio of its real size to its "
        "expected size if the full product was computed. This ratio "
        "immediately corresponds to the score, hence preferring merges "
        "where more pruning of unreachable/irrelevant states and/or perfect "
        "shrinking with bisimulation happened.");
    parser.add_option<int>("max_states", "shrink strategy option", "50000");

    options::Options options = parser.parse();
    if (parser.dry_run())
        return nullptr;
    else
        return make_shared<MergeScoringFunctionMIASM>(options);
}

static options::PluginShared<MergeScoringFunction> _plugin("miasm", _parse);
}
