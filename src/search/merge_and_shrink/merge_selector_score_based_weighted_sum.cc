#include "merge_selector_score_based_weighted_sum.h"

#include "factored_transition_system.h"
#include "types.h"

#include "../options/option_parser.h"
#include "../options/options.h"
#include "../options/plugin.h"

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
    int num_functions = merge_scoring_functions.size();
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
//    cout << "raw scores " << scores << endl;

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
//    cout << "normalized scores " << scores << endl;

    // Compute the weighted sum over all (normalized) scores.
    vector<double> weighted_sums;
    weighted_sums.reserve(num_merge_candidates);
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
        weighted_sums.push_back(weighted_sum);
    }
//    cout << "weighted sums: " << weighted_sums << endl;

    // Find the lowest scoring merge candidate.
    int best_merge_index = -1;
    double best_score = INF;
    for (int merge_index = 0; merge_index < num_merge_candidates;
         ++merge_index) {
        if (weighted_sums[merge_index] < best_score) {
            best_score = weighted_sums[merge_index];
            best_merge_index = merge_index;
        }
    }
    if (best_score == INF) {
        cerr << "All weighted summed scores are positive infinity. This "
                "cannot happen when including the total order scoring "
                "function for tie-breaking." << endl;
        utils::exit_with(utils::ExitCode::SEARCH_CRITICAL_ERROR);
    }

    int best_score_counter = 0;
    for (int merge_index = 0; merge_index < num_merge_candidates;
         ++merge_index) {
        if (weighted_sums[merge_index] == best_score) {
            ++best_score_counter;
        }
        if (best_score_counter > 1) {
            cerr << "At least two merge candidates have the exact same best "
                    "score! Did you forget to include a randomizing scoring "
                    "function for tie-breaking?" << endl;
            utils::exit_with(utils::ExitCode::SEARCH_CRITICAL_ERROR);
        }
    }

    return merge_candidates[best_merge_index];
}

string MergeSelectorScoreBasedWeightedSum::name() const {
    return "score based weighted sum";
}

void MergeSelectorScoreBasedWeightedSum::dump_specific_options() const {
    for (size_t i = 0; i < merge_scoring_functions.size(); ++i) {
        const shared_ptr<MergeScoringFunction> &scoring_function =
         merge_scoring_functions[i];
        scoring_function->dump_options();
        utils::g_log << "associated weight: " << weights[i] << endl;
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
