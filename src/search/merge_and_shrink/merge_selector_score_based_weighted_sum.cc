#include "merge_selector_score_based_weighted_sum.h"

#include "factored_transition_system.h"
#include "types.h"

#include "../options/option_parser.h"
#include "../options/options.h"
#include "../options/plugin.h"

using namespace std;

namespace merge_and_shrink {
const int MINUSINF = numeric_limits<int>::min();

double normalize_value(int min_int_score, int max_int_score, int score) {
    /* Note: we need to convert to double here because INF - MINUSINF does
       not fit into an int. */
    double min_score = min_int_score;
    double max_score = max_int_score;
    if (max_score - min_score == 0) {
        // all three scores are the same
        assert(min_score == score);
        assert(max_score == score);
        return 0;
    }
    double result = (static_cast<double>(score) - min_score) /
        (max_score - min_score);
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
}

void MergeSelectorScoreBasedWeightedSum::initialize(
    shared_ptr<AbstractTask> task) {
    for (shared_ptr<MergeScoringFunction> &scoring_function
         : merge_scoring_functions) {
        scoring_function->initialize(task);
    }
}

pair<int, int> MergeSelectorScoreBasedWeightedSum::select_merge(
    FactoredTransitionSystem &fts) const {
    // Determine the set of potential merges
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
    int num_functions = merge_scoring_functions.size();
    int num_merge_candidates = merge_candidates.size();

    // Compute all raw score values
    vector<vector<int>> unnormalized_scores(num_functions);
    for (int function_index = 0; function_index < num_functions;
         ++function_index) {
        const shared_ptr<MergeScoringFunction> &scoring_function =
            merge_scoring_functions[function_index];
        unnormalized_scores[function_index] = scoring_function->compute_scores(
            fts, merge_candidates);
    }

    // If desirge, normalize all scores.
    vector<vector<double>> normalized_scores;
    if (normalize) {
        normalized_scores.resize(num_functions);
        for (int function_index = 0; function_index < (num_functions);
             ++function_index) {
            const vector<int> &function_unn_scores =
                unnormalized_scores[function_index];
            int min_score = INF;
            int max_score = MINUSINF;
            for (int score : function_unn_scores) {
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
            for (int score: function_unn_scores) {
                double normalized_score = normalize_value(
                    min_score, max_score, score);
                function_normalized_scores.push_back(normalized_score);
            }
        }
        vector<vector<int>>().swap(unnormalized_scores);
    }

    // Compute the weighted sum over all (normalized) scores.
    vector<double> weighted_sum(num_merge_candidates, 0);
    for (int function_index = 0; function_index < num_functions;
         ++function_index) {
        int function_weight = weights[function_index];
        for (int merge_index = 0; merge_index < num_merge_candidates;
             ++merge_index) {
            weighted_sum[merge_index] += function_weight *
                (normalize ?
                 normalized_scores[function_index][merge_index] :
                 unnormalized_scores[function_index][merge_index]);
        }
    }

    // Find the loweste scoring merge candidate.
    int best_merge_index = -1;
    double best_score = INF;
    for (int merge_index = 0; merge_index < num_merge_candidates;
         ++merge_index) {
        if (weighted_sum[merge_index] < best_score) {
            best_score = weighted_sum[merge_index];
            best_merge_index = merge_index;
        }
    }

    int best_score_counter = 0;
    for (int merge_index = 0; merge_index < num_merge_candidates;
         ++merge_index) {
        if (weighted_sum[merge_index] == best_score) {
            ++best_score_counter;
        }
        if (best_score_counter > 1) {
            cerr << "At least two merge candidates have the exact same best "
                    "score! Did you forget to include a randomizing scoring "
                    "function for tie-breaking?" << endl;
            utils::exit_with(utils::ExitCode::CRITICAL_ERROR);
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
        cout << "associated weight: " << weights[i] << endl;
    }
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
        "Scores as computing by the score functions will be normalized "
        "to bin in the interval [0, 1] iff true.",
        "false");

    options::Options opts = parser.parse();

    vector<shared_ptr<MergeScoringFunction>> functions =
        opts.get_list<shared_ptr<MergeScoringFunction>>("scoring_functions");
    if (opts.contains("weights")) {
        vector<int> weights = opts.get_list<int>("weights");
        if (weights.size() != functions.size()) {
            cerr << "Number of weights differs from number of scoring "
                    "functions" << endl;
            utils::exit_with(utils::ExitCode::INPUT_ERROR);
        }
    }

    if (parser.dry_run())
        return nullptr;
    else
        return make_shared<MergeSelectorScoreBasedWeightedSum>(opts);
}

static options::PluginShared<MergeSelector> _plugin(
    "score_based_weighted_sum", _parse);
}