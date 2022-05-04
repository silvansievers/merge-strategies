#include "merge_selector_score_based_weighted_sum.h"

#include "factored_transition_system.h"
#include "types.h"

#include "../options/option_parser.h"
#include "../options/options.h"
#include "../options/plugin.h"

#include "../utils/hash.h"
#include "../utils/logging.h"

using namespace std;

namespace merge_and_shrink {
double normalize_value(double min_score, double max_score, double score) {
    if (max_score - min_score == 0) {
        // all three scores are the same
        assert(min_score == score);
        assert(max_score == score);
        return 0;
    }
    double result = (score - min_score) / (max_score - min_score);
    assert(result >= 0);
    assert(result <= 1);
    return result;
}

MergeSelectorScoreBasedWeightedSum::MergeSelectorScoreBasedWeightedSum(
    const options::Options &options)
    : merge_scoring_functions(
          options.get_list<shared_ptr<MergeScoringFunction>>(
              "scoring_functions")),
      normalize(options.get<bool>("normalize")) {
    if (options.contains("weights")) {
        weights = options.get_list<int>("weights");
    } else {
        weights.resize(merge_scoring_functions.size(), 1);
    }
    if (normalize) {
        // TODO: normalizing infinite values doesn't work.
        cerr << "Normalizing scores currently doesn't work." << endl;
        utils::exit_with(utils::ExitCode::SEARCH_CRITICAL_ERROR);
    }
    if (merge_scoring_functions.back()->get_name() != "total order") {
        cerr << "This merge selector requires using the total order scoring "
            "function as the last in the list for tie-breaking in cases "
            "where two or more merge candidates have the same best "
            "weighted sum as score." << endl;
        utils::exit_with(utils::ExitCode::SEARCH_INPUT_ERROR);
    }
}

void MergeSelectorScoreBasedWeightedSum::initialize(
    const TaskProxy &task_proxy) {
    for (shared_ptr<MergeScoringFunction> &scoring_function
         : merge_scoring_functions) {
        scoring_function->initialize(task_proxy);
    }
}

pair<int, int> MergeSelectorScoreBasedWeightedSum::select_merge(
    const FactoredTransitionSystem &fts,
    const vector<int> &indices_subset) const {
    vector<pair<int, int>> merge_candidates =
        compute_merge_candidates(fts, indices_subset);
    // The last scoring function has to be "total order" for tie-breaking.
    int num_functions = merge_scoring_functions.size() - 1;
    int num_merge_candidates = merge_candidates.size();

    // Compute all raw score values
    vector<vector<double>> scores(num_functions);
    for (int function_index = 0; function_index < num_functions;
         ++function_index) {
        const shared_ptr<MergeScoringFunction> &scoring_function =
            merge_scoring_functions[function_index];
        scores[function_index] = scoring_function->compute_scores(
            fts, merge_candidates);
    }
//    cout << "raw scores" << endl;
//    for (int cand_ind = 0; cand_ind < num_merge_candidates; ++cand_ind) {
//        cout << "candidate " << merge_candidates[cand_ind].first
//             << ", " << merge_candidates[cand_ind].second << ": ";
//        for (int function = 0; function < num_functions; ++function) {
//             cout << scores[function][cand_ind] << ", ";
//        }
//        cout << endl;
//    }

    // If desired, normalize all scores.
    if (normalize) {
        vector<vector<double>> normalized_scores;
        normalized_scores.resize(num_functions);
        for (int function_index = 0; function_index < (num_functions);
             ++function_index) {
            const vector<double> &raw_scores = scores[function_index];
            double min_score = INF;
            double max_score = MINUSINF;
            for (double score : raw_scores) {
                if (score < min_score) {
                    min_score = score;
                }
                if (score > max_score) {
                    max_score = score;
                }
            }

            vector<double> &function_normalized_scores =
                normalized_scores[function_index];
            function_normalized_scores.reserve(num_merge_candidates);
            for (double score: raw_scores) {
                double normalized_score = normalize_value(
                    min_score, max_score, score);
                function_normalized_scores.push_back(normalized_score);
            }
        }
        scores.swap(normalized_scores);
    }
//    cout << "normalized scores" << endl;
//    for (int cand_ind = 0; cand_ind < num_merge_candidates; ++cand_ind) {
//        cout << "candidate " << merge_candidates[cand_ind].first
//             << ", " << merge_candidates[cand_ind].second << ": ";
//        for (int function = 0; function < num_functions; ++function) {
//            cout << scores[function][cand_ind] << ", ";
//        }
//        cout << endl;
//    }

    // Compute the weighted sum over all (normalized) scores.
    vector<double> weighted_summed_scores;
    weighted_summed_scores.reserve(num_merge_candidates);
    for (int merge_index = 0; merge_index < num_merge_candidates;
         ++merge_index) {
        double weighted_sum = 0;
        // First check if any weight is MINUSINF which takes precedence.
        for (int function_index = 0; function_index < num_functions;
             ++function_index) {
            double score = scores[function_index][merge_index];
            if (score == MINUSINF) {
                weighted_sum = MINUSINF;
                break;
            }
        }
        if (weighted_sum == 0) {
            bool all_scores_infinity = true;
            for (int function_index = 0; function_index < num_functions;
                 ++function_index) {
                double score = scores[function_index][merge_index];
                int function_weight = weights[function_index];
                if (score != INF) {
                    weighted_sum += function_weight * score;
                    all_scores_infinity = false;
                }
            }
            if (all_scores_infinity) {
                weighted_sum = INF;
            }
        }
        weighted_summed_scores.push_back(weighted_sum);
    }
//    cout << "weighted summed scores" << endl;
//    for (int cand_ind = 0; cand_ind < num_merge_candidates; ++cand_ind) {
//        cout << "candidate " << merge_candidates[cand_ind].first
//             << ", " << merge_candidates[cand_ind].second << ": "
//             << weighted_summed_scores[cand_ind] << endl;
//    }

    // Find the lowest scoring merge candidate.
    double best_score = INF;
    vector<int> best_merge_indices;
    for (int merge_index = 0; merge_index < num_merge_candidates;
         ++merge_index) {
        int score = weighted_summed_scores[merge_index];
        if (score < best_score) {
            best_score = score;
            best_merge_indices.clear();
            best_merge_indices.push_back(merge_index);
        } else if (score == best_score) {
            best_merge_indices.push_back(merge_index);
        }
    }
    if (best_score == INF) {
        cerr << "All weighted summed scores are positive infinity. This "
            "cannot happen when including the total order scoring "
            "function for tie-breaking." << endl;
        utils::exit_with(utils::ExitCode::SEARCH_CRITICAL_ERROR);
    }

    int best_merge_index = -1;
    if (best_merge_indices.size() == 1) {
        best_merge_index = best_merge_indices.front();
    } else {
        const shared_ptr<MergeScoringFunction> &total_order_scoring_function =
            merge_scoring_functions.back();
        vector<double> tie_breaking_scores =
            total_order_scoring_function->compute_scores(fts, merge_candidates);
        double best_score = INF;
        for (int merge_index : best_merge_indices) {
            int score = tie_breaking_scores[merge_index];
            assert(score != best_score);
            if (score < best_score) {
                best_score = score;
                best_merge_index = merge_index;
            }
        }
    }

    return merge_candidates[best_merge_index];
}

string MergeSelectorScoreBasedWeightedSum::name() const {
    return "score based weighted sum";
}

void MergeSelectorScoreBasedWeightedSum::dump_selector_specific_options(
    utils::LogProxy &log) const {
    for (size_t i = 0; i < merge_scoring_functions.size(); ++i) {
        const shared_ptr<MergeScoringFunction> &scoring_function =
            merge_scoring_functions[i];
        scoring_function->dump_options(log);
        log << "associated weight: " << weights[i] << endl;
    }
}

bool MergeSelectorScoreBasedWeightedSum::requires_init_distances() const {
    for (const shared_ptr<MergeScoringFunction> &scoring_function
         : merge_scoring_functions) {
        if (scoring_function->requires_init_distances()) {
            return true;
        }
    }
    return false;
}

bool MergeSelectorScoreBasedWeightedSum::requires_goal_distances() const {
    for (const shared_ptr<MergeScoringFunction> &scoring_function
         : merge_scoring_functions) {
        if (scoring_function->requires_goal_distances()) {
            return true;
        }
    }
    return false;
}

static shared_ptr<MergeSelector>_parse(options::OptionParser &parser) {
    parser.document_synopsis(
        "Score based weighted sum merge selector",
        "This merge selector has a list of scoring functions. The final score "
        "assigned to each merge candidate is the weighted sum over all "
        "scores of all scoring functions, normalized to be in the interval "
        "[0, 1]. The merge candidate with the minimum score is selected.");
    parser.add_list_option<shared_ptr<MergeScoringFunction>>(
        "scoring_functions",
        "The list of scoring functions used to compute scores for candidates.");
    parser.add_list_option<int>(
        "weights",
        "The list of weights to be used for the scoring functions. The list"
        "must either be empty, in which case the same weight (1) for each"
        "scoring function is assumed, or its size must match the number of"
        "given scoring functions, assigning a weight to each.",
        options::OptionParser::NONE);
    parser.add_option<bool>(
        "normalize",
        "Scores as compute by the score functions will be normalized "
        "to be in the interval [0, 1] iff true.",
        "false");

    options::Options opts = parser.parse();
    if (parser.help_mode()) {
        return nullptr;
    }

    vector<shared_ptr<MergeScoringFunction>> functions =
        opts.get_list<shared_ptr<MergeScoringFunction>>("scoring_functions");
    if (opts.contains("weights")) {
        vector<int> weights = opts.get_list<int>("weights");
        if (weights.size() != functions.size()) {
            cerr << "Number of weights differs from number of scoring "
                "functions" << endl;
            utils::exit_with(utils::ExitCode::SEARCH_INPUT_ERROR);
        }
    }

    if (parser.dry_run())
        return nullptr;
    else
        return make_shared<MergeSelectorScoreBasedWeightedSum>(opts);
}

static options::Plugin<MergeSelector> _plugin(
    "score_based_weighted_sum", _parse);
}
