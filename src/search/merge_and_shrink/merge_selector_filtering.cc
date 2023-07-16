#include "merge_selector_filtering.h"

#include "factored_transition_system.h"
#include "merge_scoring_function.h"

#include "../task_proxy.h"

#include "../plugins/plugin.h"

#include <cassert>

using namespace std;

namespace merge_and_shrink {
MergeSelectorFiltering::MergeSelectorFiltering(
    const plugins::Options &options)
    : merge_scoring_functions(
          options.get_list<shared_ptr<MergeScoringFunction>>(
              "scoring_functions")),
      iterations_with_tiebreaking(0),
      total_tiebreaking_pair_count(0),
      num_candidates(0) {
}

shared_ptr<MergeCandidate> MergeSelectorFiltering::get_candidate(
    int index1, int index2) const {
    assert(utils::in_bounds(index1, merge_candidates_by_indices));
    assert(utils::in_bounds(index2, merge_candidates_by_indices[index1]));
    if (merge_candidates_by_indices[index1][index2] == nullptr) {
        merge_candidates_by_indices[index1][index2] =
            make_shared<MergeCandidate>(num_candidates, index1, index2);
        ++num_candidates;
    }
    return merge_candidates_by_indices[index1][index2];
}

static vector<shared_ptr<MergeCandidate>> get_remaining_candidates(
    const vector<shared_ptr<MergeCandidate>> &merge_candidates,
    const vector<double> &scores) {
    assert(merge_candidates.size() == scores.size());
    double best_score = INF;
    for (double score : scores) {
        if (score < best_score) {
            best_score = score;
        }
    }

    vector<shared_ptr<MergeCandidate>> result;
    for (size_t i = 0; i < scores.size(); ++i) {
        if (scores[i] == best_score) {
            result.push_back(merge_candidates[i]);
        }
    }
    return result;
}

pair<int, int> MergeSelectorFiltering::select_merge(
    const FactoredTransitionSystem &fts,
    const vector<int> &indices_subset) const {
    vector<shared_ptr<MergeCandidate>> merge_candidates;
    if (indices_subset.empty()) {
        for (int ts_index1 = 0; ts_index1 < fts.get_size(); ++ts_index1) {
            if (fts.is_active(ts_index1)) {
                for (int ts_index2 = ts_index1 + 1; ts_index2 < fts.get_size();
                     ++ts_index2) {
                    if (fts.is_active(ts_index2)) {
                        merge_candidates.push_back(get_candidate(ts_index1, ts_index2));
                    }
                }
            }
        }
    } else {
        assert(indices_subset.size() > 1);
        for (size_t i = 0; i < indices_subset.size(); ++i) {
            int ts_index1 = indices_subset[i];
            assert(fts.is_active(ts_index1));
            for (size_t j = i + 1; j < indices_subset.size(); ++j) {
                int ts_index2 = indices_subset[j];
                assert(fts.is_active(ts_index2));
                merge_candidates.push_back(get_candidate(ts_index1, ts_index2));
            }
        }
    }

    for (size_t i = 0; i < merge_scoring_functions.size(); ++i) {
        const shared_ptr<MergeScoringFunction> &scoring_function =
            merge_scoring_functions[i];
        vector<double> scores = scoring_function->compute_scores_caching(
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

    return make_pair(merge_candidates.front()->index1, merge_candidates.front()->index2);
}

void MergeSelectorFiltering::initialize(const TaskProxy &task_proxy) {
    for (shared_ptr<MergeScoringFunction> &scoring_function
         : merge_scoring_functions) {
        scoring_function->initialize(task_proxy);
    }
    int num_variables = task_proxy.get_variables().size();
    int max_factor_index = 2 * num_variables - 1;

    merge_candidates_by_indices.resize(
        max_factor_index,
        vector<shared_ptr<MergeCandidate>>(max_factor_index, nullptr));
}

string MergeSelectorFiltering::name() const {
    return "score based filtering";
}

void MergeSelectorFiltering::dump_selector_specific_options(
    utils::LogProxy &log) const {
    if (log.is_at_least_normal()) {
        for (const shared_ptr<MergeScoringFunction> &scoring_function
             : merge_scoring_functions) {
            scoring_function->dump_options(log);
        }
    }
}

bool MergeSelectorFiltering::requires_init_distances() const {
    for (const shared_ptr<MergeScoringFunction> &scoring_function
         : merge_scoring_functions) {
        if (scoring_function->requires_init_distances()) {
            return true;
        }
    }
    return false;
}

bool MergeSelectorFiltering::requires_goal_distances() const {
    for (const shared_ptr<MergeScoringFunction> &scoring_function
         : merge_scoring_functions) {
        if (scoring_function->requires_goal_distances()) {
            return true;
        }
    }
    return false;
}

class MergeSelectorFilteringFeature : public plugins::TypedFeature<MergeSelector, MergeSelectorFiltering> {
public:
    MergeSelectorFilteringFeature() : TypedFeature("filtering") {
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

static plugins::FeaturePlugin<MergeSelectorFilteringFeature> _plugin;
}
