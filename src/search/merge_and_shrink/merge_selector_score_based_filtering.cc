#include "merge_selector_score_based_filtering.h"

#include "factored_transition_system.h"
#include "merge_scoring_function.h"

#include "../plugins/plugin.h"

#include <cassert>

using namespace std;

namespace merge_and_shrink {
MergeSelectorScoreBasedFiltering::MergeSelectorScoreBasedFiltering(
    const plugins::Options &options)
    : merge_scoring_functions(
          options.get_list<shared_ptr<MergeScoringFunction>>(
              "scoring_functions")),
      iterations_with_tiebreaking(0),
      total_tiebreaking_pair_count(0) {
}

vector<pair<int, int>> MergeSelectorScoreBasedFiltering::get_remaining_candidates(
    const vector<pair<int, int>> &merge_candidates,
    const vector<double> &scores) const {
    assert(merge_candidates.size() == scores.size());
    double best_score = INF;
    for (double score : scores) {
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
    const FactoredTransitionSystem &fts,
    const vector<int> &indices_subset) const {
    vector<pair<int, int>> merge_candidates =
        compute_merge_candidates(fts, indices_subset);

    for (size_t i = 0; i < merge_scoring_functions.size(); ++i) {
        const shared_ptr<MergeScoringFunction> &scoring_function =
            merge_scoring_functions[i];
        vector<double> scores = scoring_function->compute_scores(
            fts, merge_candidates);
        if ((scoring_function->get_name() == "total order" ||
             scoring_function->get_name() == "single random") &&
            merge_candidates.size() != 1) {
            ++iterations_with_tiebreaking;
            total_tiebreaking_pair_count += merge_candidates.size();
        }
        int old_size = merge_candidates.size();
        merge_candidates = get_remaining_candidates(merge_candidates, scores);
        int new_size = merge_candidates.size();
        if (new_size - old_size == 0 &&
            i != merge_scoring_functions.size() - 1 &&
            scoring_function->get_name() == "goal relevance" &&
            merge_scoring_functions[i + 1]->get_name() == "dfp" &&
            scores.front() == INF) {
            /*
              "skip" computation of dfp if no goal relevance pair has been
              found (which is the case if no candidates have been filtered,
              because all scores are identical, and the scores must be
              infinity) to mimic previous MergeDFP behavior.
            */
            ++i;
            continue;
        }
        if (merge_candidates.size() == 1) {
            break;
        }
    }

    if (merge_candidates.size() > 1) {
        cerr << "More than one merge candidate remained after computing all "
            "scores! Did you forget to include a uniquely tie-breaking "
            "scoring function, e.g. total_order or single_random?" << endl;
        utils::exit_with(utils::ExitCode::SEARCH_CRITICAL_ERROR);
    }

    return merge_candidates.front();
}

void MergeSelectorScoreBasedFiltering::initialize(const TaskProxy &task_proxy) {
    for (shared_ptr<MergeScoringFunction> &scoring_function
         : merge_scoring_functions) {
        scoring_function->initialize(task_proxy);
    }
}

string MergeSelectorScoreBasedFiltering::name() const {
    return "score based filtering";
}

void MergeSelectorScoreBasedFiltering::dump_selector_specific_options(
    utils::LogProxy &log) const {
    if (log.is_at_least_normal()) {
        for (const shared_ptr<MergeScoringFunction> &scoring_function
             : merge_scoring_functions) {
            scoring_function->dump_options(log);
        }
    }
}

bool MergeSelectorScoreBasedFiltering::requires_init_distances() const {
    for (const shared_ptr<MergeScoringFunction> &scoring_function
         : merge_scoring_functions) {
        if (scoring_function->requires_init_distances()) {
            return true;
        }
    }
    return false;
}

bool MergeSelectorScoreBasedFiltering::requires_goal_distances() const {
    for (const shared_ptr<MergeScoringFunction> &scoring_function
         : merge_scoring_functions) {
        if (scoring_function->requires_goal_distances()) {
            return true;
        }
    }
    return false;
}

class MergeSelectorScoreBasedFilteringFeature : public plugins::TypedFeature<MergeSelector, MergeSelectorScoreBasedFiltering> {
public:
    MergeSelectorScoreBasedFilteringFeature() : TypedFeature("score_based_filtering") {
        document_title("Score based filtering merge selector");
        document_synopsis(
            "This merge selector has a list of scoring functions, which are used "
            "iteratively to compute scores for merge candidates, keeping the best "
            "ones (with minimal scores) until only one is left.");

        add_list_option<shared_ptr<MergeScoringFunction>>(
            "scoring_functions",
            "The list of scoring functions used to compute scores for candidates.");
    }
};

static plugins::FeaturePlugin<MergeSelectorScoreBasedFilteringFeature> _plugin;
}
