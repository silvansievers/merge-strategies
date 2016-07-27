#include "merge_selector_score_based_filtering.h"

#include "factored_transition_system.h"
#include "merge_scoring_function.h"

#include "../options/option_parser.h"
#include "../options/options.h"
#include "../options/plugin.h"

#include <cassert>

using namespace std;

namespace merge_and_shrink {
MergeSelectorScoreBasedFiltering::MergeSelectorScoreBasedFiltering(
    const options::Options &options)
    : merge_scoring_functions(
        options.get_list<shared_ptr<MergeScoringFunction>>(
            "scoring_functions")),
      iterations_with_tiebreaking(0),
      total_tiebreaking_pair_count(0) {
}

MergeSelectorScoreBasedFiltering::MergeSelectorScoreBasedFiltering(
    vector<shared_ptr<MergeScoringFunction>> scoring_functions)
    : merge_scoring_functions(move(scoring_functions)),
      iterations_with_tiebreaking(0),
      total_tiebreaking_pair_count(0) {
}

vector<pair<int, int>> MergeSelectorScoreBasedFiltering::get_remaining_candidates(
    const vector<pair<int, int>> &merge_candidates,
    const vector<int> &scores) const {
    assert(merge_candidates.size() == scores.size());
    int best_score = INF;
    for (int score : scores) {
        if (score < best_score) {
            best_score = score;
        }
    }

    vector<pair<int, int>> result;
    for (size_t i = 0; i < scores.size(); ++i) {
        if (scores[i] == best_score) {
            result.push_back(merge_candidates[i]);
        }
    }
    return result;
}

pair<int, int> MergeSelectorScoreBasedFiltering::select_merge(
    FactoredTransitionSystem &fts) const {
    vector<pair<int, int>> merge_candidates;
    for (int ts_index1 = 0; ts_index1 < fts.get_size(); ++ts_index1) {
        if (fts.is_active(ts_index1)) {
            for (int ts_index2 = ts_index1 + 1; ts_index2 < fts.get_size();
                 ++ts_index2) {
                if (fts.is_active(ts_index2)) {
                    merge_candidates.emplace_back(ts_index1, ts_index2);
                }
            }
        }
    }

    for (const shared_ptr<MergeScoringFunction> &scoring_function :
         merge_scoring_functions) {
        vector<int> scores = scoring_function->compute_scores(fts,
                                                              merge_candidates);
        merge_candidates = get_remaining_candidates(merge_candidates, scores);

        if (scoring_function->name() == "dfp" && merge_candidates.size() != 1) {
            ++iterations_with_tiebreaking;
            total_tiebreaking_pair_count += merge_candidates.size();
        }

        if (merge_candidates.size() == 1) {
            break;
        }
    }

    if (merge_candidates.size() > 1) {
        cerr << "More than one merge candidate remained after computing all "
                "scores! Did you forget to include a randomizing scoring "
                "function for tie-breaking?" << endl;
        utils::exit_with(utils::ExitCode::CRITICAL_ERROR);
    }

    return merge_candidates.front();
}

// TODO: code duplication!
pair<int, int> MergeSelectorScoreBasedFiltering::select_merge_dfp_sccs(
    FactoredTransitionSystem &fts,
    const vector<int> &indices_subset) {
    vector<pair<int, int>> merge_candidates;
    for (size_t i = 0; i < indices_subset.size(); ++i) {
        int ts_index1 = indices_subset[i];
        assert(fts.is_active(ts_index1));
        for (size_t j = i + 1; j < indices_subset.size(); ++j) {
            int ts_index2 = indices_subset[j];
            assert(fts.is_active(ts_index2));
            merge_candidates.emplace_back(ts_index1, ts_index2);
        }
    }

    for (size_t i = 0; i < merge_scoring_functions.size(); ++i) {
        shared_ptr<MergeScoringFunction> &scoring_function =
            merge_scoring_functions[i];
        vector<int> scores = scoring_function->compute_scores(fts,
                                                              merge_candidates);
        int old_size = merge_candidates.size();
        merge_candidates = get_remaining_candidates(merge_candidates, scores);
        int new_size = merge_candidates.size();
        if (new_size - old_size == 0 &&
            scoring_function->name() == "goal relevance" &&
            scores.front() == INF) {
            /*
              "skip" computation of dfp if no goal relevance pair has been
              found (which is the case if no candidates have been filtered,
              because all scores are identical, and the scores must be
              infinity) to mimic previous MergeDFP behavior.
            */
            assert(i < merge_scoring_functions.size() - 1);
            assert(merge_scoring_functions[i + 1]->name() == "dfp");
            ++i;
//            cout << "skipping computation of dfp scoring function" << endl;
            continue;
        }

        if (scoring_function->name() == "dfp" && merge_candidates.size() != 1) {
            ++iterations_with_tiebreaking;
            total_tiebreaking_pair_count += merge_candidates.size();
        }

        if (merge_candidates.size() == 1) {
            break;
        }
    }

    if (merge_candidates.size() > 1) {
        cerr << "More than one merge candidate remained after computing all "
                "scores! Did you forget to include a randomizing scoring "
                "function for tie-breaking?" << endl;
        utils::exit_with(utils::ExitCode::CRITICAL_ERROR);
    }

    return merge_candidates.front();
}

void MergeSelectorScoreBasedFiltering::initialize(shared_ptr<AbstractTask> task) {
    for (shared_ptr<MergeScoringFunction> &scoring_function
         : merge_scoring_functions) {
        scoring_function->initialize(task);
    }
}

string MergeSelectorScoreBasedFiltering::name() const {
    return "score based filtering";
}

void MergeSelectorScoreBasedFiltering::dump_specific_options() const {
    for (const shared_ptr<MergeScoringFunction> &scoring_function
         : merge_scoring_functions) {
        scoring_function->dump_options();
    }
}

static shared_ptr<MergeSelector>_parse(options::OptionParser &parser) {
    parser.document_synopsis(
        "Score based filtering merge selector",
        "This merge selector has a list of scoring functions, which are used "
        "iteratively to compute scores for merge candidates, keeping the best "
        "ones until only one is left.");
    parser.add_list_option<shared_ptr<MergeScoringFunction>>(
        "scoring_functions",
        "The list of scoring functions used to compute scores for candidates.");

    options::Options opts = parser.parse();
    if (parser.dry_run())
        return nullptr;
    else
        return make_shared<MergeSelectorScoreBasedFiltering>(opts);
}

static options::PluginShared<MergeSelector> _plugin("score_based_filtering", _parse);
}
