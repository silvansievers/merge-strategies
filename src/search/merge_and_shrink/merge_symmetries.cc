#include "merge_symmetries.h"

#include "factored_transition_system.h"
#include "merge_selector.h"
#include "merge_tree.h"
#include "transition_system.h"

#include "symmetries/symmetry_group.h"

#include "../options/options.h"

#include <algorithm>
#include <iomanip>

using namespace std;

namespace merge_and_shrink {
MergeSymmetries::MergeSymmetries(
    FactoredTransitionSystem &fts,
    const options::Options &options,
    int num_merges,
    unique_ptr<MergeTree> merge_tree,
    shared_ptr<MergeSelector> merge_selector)
    : MergeStrategy(fts),
      options(options),
      num_merges(num_merges),
      merge_tree(move(merge_tree)),
      merge_selector(merge_selector),
      max_bliss_iterations(options.get<int>("max_bliss_iterations")),
      bliss_call_time_limit(options.get<int>("bliss_call_time_limit")),
      bliss_remaining_time_budget(options.get<int>("bliss_total_time_budget")),
      iteration_counter(0),
      number_of_applied_symmetries(0),
      bliss_limit_reached(false),
      pure_fallback_strategy(true) {
}

void MergeSymmetries::dump_statistics() {
    if (!bliss_times.empty()) {
        sort(bliss_times.begin(), bliss_times.end());
        double summed_up_bliss_times = 0;
        size_t total_bliss_calls = bliss_times.size();
        cout << "Total bliss calls: " << total_bliss_calls << endl;
        for (size_t i = 0; i < total_bliss_calls; ++i) {
            summed_up_bliss_times += bliss_times[i];
        }
        double average_bliss_time = summed_up_bliss_times
                / (double) total_bliss_calls;
        cout << setprecision(5) << fixed << "Average bliss time: " << average_bliss_time << endl;
        double median_bliss_time;
        if (total_bliss_calls % 2 == 0) {
            size_t lower_median_index = (total_bliss_calls - 1) / 2;
            size_t higher_median_index = lower_median_index + 1;
            median_bliss_time = (bliss_times[lower_median_index]
                                 + bliss_times[higher_median_index])
                    / 2.0;
        } else {
            size_t median_index = total_bliss_calls / 2;
            median_bliss_time = bliss_times[median_index];
        }
        cout << setprecision(5) << fixed << "Median bliss time: " << median_bliss_time << endl;
    }
}

pair<int, int> MergeSymmetries::get_next() {
    ++iteration_counter;

    if (!bliss_limit_reached && iteration_counter <= max_bliss_iterations && merge_order.empty()) {
        double time_limit = 0.0;
        if (bliss_remaining_time_budget) {
            // we round everything down to the full second
            time_limit = bliss_remaining_time_budget;
        } else if (bliss_call_time_limit) {
            time_limit = bliss_call_time_limit;
        }
        cout << "Setting bliss time limit to " << time_limit << endl;
        options.set<double>("bliss_time_limit", time_limit);
        SymmetryGroup symmetry_group(options);
        bool applied_symmetries =
            symmetry_group.find_and_apply_symmetries(fts, merge_order);
        if (applied_symmetries) {
            ++number_of_applied_symmetries;
        }
        if (symmetry_group.is_bliss_limit_reached()) {
            bliss_limit_reached = true;
        }
        if (pure_fallback_strategy && (!merge_order.empty())) {
            pure_fallback_strategy = false;
            cout << "not pure fallback strategy anymore" << endl;
        }
        double bliss_time = symmetry_group.get_bliss_time();
        bliss_times.push_back(bliss_time);
        if (bliss_remaining_time_budget) {
            bliss_remaining_time_budget -= bliss_time;
            if (bliss_remaining_time_budget <= 0) {
                assert(bliss_limit_reached);
                // in case of measurement inaccuracies, set bliss limit
                // reached to true in any case.
                bliss_limit_reached = true;
            }
            cout << "Remaining bliss time budget " << bliss_remaining_time_budget << endl;
        }
        cout << "Number of applied symmetries: " << number_of_applied_symmetries << endl;
    }

    if (iteration_counter == num_merges) { // TODO: check correctness
        dump_statistics();
    }

    int next_merge_index = fts.get_size();

    // No symmetries merge order -- use fallback strategy.
    if (merge_order.empty()) {
        if (merge_selector) {
            return merge_selector->select_merge(fts);
        } else {
            assert(merge_tree);
            return merge_tree->get_next_merge(next_merge_index);
        }
    }

    pair<int, int> next_merge = merge_order.front();
    if (!fts.is_active(next_merge.first) || !fts.is_active(next_merge.second)) {
        cerr << "Problem with the merge strategy: invalid indices" << endl;
        utils::exit_with(utils::ExitCode::CRITICAL_ERROR);
    }
    merge_order.erase(merge_order.begin());

    // Update the merge tree, if present.
    if (merge_tree) {
        merge_tree->update(
            next_merge,
            next_merge_index);
    }

    return next_merge;
}
}
