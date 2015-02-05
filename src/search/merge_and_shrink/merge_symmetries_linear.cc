#include "merge_symmetries_linear.h"

#include "symmetries/symmetries.h"

#include "../globals.h"
#include "../plugin.h"
#include "../variable_order_finder.h"

#include <algorithm>
#include <iomanip>

using namespace std;

MergeSymmetriesLinear::MergeSymmetriesLinear(const Options &options_)
    : MergeStrategy(),
      options(options_),
      max_bliss_iterations(options.get<int>("max_bliss_iterations")),
      bliss_call_time_limit(options.get<int>("bliss_call_time_limit")),
      bliss_remaining_time_budget(options.get<int>("bliss_total_time_budget")),
      iteration_counter(0),
      number_of_applied_symmetries(0),
      bliss_limit_reached(false),
      only_applied_dfp(true) {
    VariableOrderFinder order(VariableOrderType(options.get_enum("variable_order")));
    linear_merge_order.reserve(g_variable_domain.size());
    while (!order.done()) {
        linear_merge_order.push_back(order.next());
    }
    cout << "linear merge order: " << linear_merge_order << endl;
}

void MergeSymmetriesLinear::dump_statistics() {
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

void MergeSymmetriesLinear::dump_strategy_specific_options() const {
    cout << "Options for merge symmetries:" << endl;
    cout << "    symmetries for shrinking: ";
    int symmetries_for_shrinking = options.get_enum("symmetries_for_shrinking");
    switch (symmetries_for_shrinking) {
        case 0: {
            cout << "none";
            break;
        }
        case 1: {
            cout << "atomic";
            break;
        }
        case 2: {
            cout << "local";
            break;
        }
    }
    cout << endl;
    cout << "    symmetries for merging: ";
    int symmetries_for_merging = options.get_enum("symmetries_for_merging");
    switch (symmetries_for_merging) {
        case 0: {
            cout << "none";
            break;
        }
        case 1: {
            cout << "smallest";
            break;
        }
        case 2: {
            cout << "largest";
            break;
        }
    }
    cout << endl;
    if (symmetries_for_merging) {
        cout << "    external merging: ";
        switch (options.get_enum("external_merging")) {
            case 0: {
                cout << "merge for atomic symmetry";
                break;
            }
            case 1: {
                cout << "merge for local symmetry";
                break;
            }
        }
        cout << endl;
        cout << "    internal merging: ";
        switch (options.get_enum("internal_merging")) {
            case 0: {
                cout << "linear";
                break;
            }
            case 1: {
                cout << "non linear";
                break;
            }
        }
        cout << endl;
    }
    cout << "    maxium number of m&s iterations with bliss: "
         << max_bliss_iterations << endl;
    cout << "    time limit for single bliss calls (0 means unlimited): "
         << bliss_call_time_limit << endl;
    cout << "    total time budget for bliss (0 means unlimited): "
         << options.get<int>("bliss_total_time_budget") << endl;
    cout << "    stop searching for symmetries once no symmetry was found: "
         << (options.get<bool>("stop_after_no_symmetries") ? "yes" : "no") << endl;
    cout << "    build stabilized pdg: "
         << (options.get<bool>("build_stabilized_pdg") ? "yes" : "no") << endl;
}

pair<int, int> MergeSymmetriesLinear::get_next(const std::vector<TransitionSystem *> &all_abstractions) {
    assert(!done());
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
        Symmetries symmetries(options);
        bool applied_symmetries =
                symmetries.find_and_apply_symmetries(all_abstractions, merge_order);
        if (applied_symmetries) {
            ++number_of_applied_symmetries;
        }
        if (symmetries.is_bliss_limit_reached()) {
            bliss_limit_reached = true;
        }
        if (only_applied_dfp && (applied_symmetries || !merge_order.empty())) {
            only_applied_dfp = false;
            cout << "not pure DFP anymore" << endl;
        }
        double bliss_time = symmetries.get_bliss_time();
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

    if (remaining_merges == 1) {
        dump_statistics();
    }

    if (merge_order.empty()) {
        // Linear Strategy
        int first = linear_merge_order[0];
        linear_merge_order.erase(linear_merge_order.begin());
        int second = linear_merge_order[0];
        linear_merge_order[0] = all_abstractions.size();
        cout << "Next pair (linear strategy): " << first << ", " << second << endl;
        --remaining_merges;
        return make_pair(first, second);
    }

    pair<int, int> next_merge = merge_order.front();
    if (!all_abstractions[next_merge.first] || !all_abstractions[next_merge.second]) {
        cerr << "Problem with the merge strategy: invalid indices" << endl;
        exit_with(EXIT_CRITICAL_ERROR);
    }
    merge_order.erase(merge_order.begin());

    bool found_entry = false;
    size_t counter = 0;
    for (vector<int>::iterator it = linear_merge_order.begin();
         it != linear_merge_order.end(); ++it) {
        if (!found_entry && (*it == next_merge.first ||
                             *it == next_merge.second)) {
            //cout << "updating entry " << counter << " (" << *it << ")to " << all_abstractions.size() << endl;
            found_entry = true;
            linear_merge_order[counter] = all_abstractions.size();
        } else if (found_entry && (*it == next_merge.first ||
                                   *it == next_merge.second)) {
            //cout << "erasing entry " << counter << " (" << *it << ")" << endl;
            linear_merge_order.erase(it);
            break;
        }
        ++counter;
    }
    assert(found_entry);

    --remaining_merges;
    return next_merge;
}

string MergeSymmetriesLinear::name() const {
    return "symmetries";
}

static MergeStrategy *_parse(OptionParser &parser) {
    vector<string> merge_strategies;
    merge_strategies.push_back("CG_GOAL_LEVEL");
    merge_strategies.push_back("CG_GOAL_RANDOM");
    merge_strategies.push_back("GOAL_CG_LEVEL");
    merge_strategies.push_back("RANDOM");
    merge_strategies.push_back("LEVEL");
    merge_strategies.push_back("REVERSE_LEVEL");
    parser.add_enum_option("variable_order", merge_strategies,
                           "the order in which atomic abstractions are merged",
                           "CG_GOAL_LEVEL");

    parser.add_option<int>("max_bliss_iterations", "maximum ms iteration until "
                           "which bliss is allowed to run.",
                           "infinity");
    parser.add_option<int>("bliss_call_time_limit", "time in seconds one bliss "
                           "run is allowed to last at most (0 means no limit)",
                           "0");
    parser.add_option<int>("bliss_total_time_budget", "time in seconds bliss is "
                           "allowed to run overall (0 means no limit)",
                           "0");
    parser.add_option<bool>("stop_after_no_symmetries", "stop calling bliss "
                            "after unsuccesfull previous bliss call.",
                           "False");
    vector<string> symmetries_for_shrinking;
    symmetries_for_shrinking.push_back("NO_SHRINKING");
    symmetries_for_shrinking.push_back("ATOMIC");
    symmetries_for_shrinking.push_back("LOCAL");
    parser.add_enum_option("symmetries_for_shrinking",
                           symmetries_for_shrinking,
                           "choose the type of symmetries used for shrinking: "
                           "no shrinking, "
                           "only atomic symmetries, "
                           "local symmetries.",
                           "NO_SHRINKING");
    vector<string> symmetries_for_merging;
    symmetries_for_merging.push_back("NO_MERGING");
    symmetries_for_merging.push_back("SMALLEST");
    symmetries_for_merging.push_back("LARGEST");
    parser.add_enum_option("symmetries_for_merging",
                           symmetries_for_merging,
                           "choose the type of symmetries that should determine "
                           "the set of abstractions to be merged: "
                           "the smallest or the largest",
                           "SMALLEST");
    vector<string> external_merging;
    external_merging.push_back("MERGE_FOR_ATOMIC");
    external_merging.push_back("MERGE_FOR_LOCAL");
    parser.add_enum_option("external_merging",
                           external_merging,
                           "choose the set of abstractions to be merged: "
                           "merge for atomic: merge all abstractions affected "
                           "by the chosen symmetry, or "
                           "merge for local: merge only the abstractions "
                           "mapped (in cycles) to others. only merge every "
                           "cycle separately.",
                           "MERGE_FOR_ATOMIC");
    vector<string> internal_merging;
    internal_merging.push_back("LINEAR");
    internal_merging.push_back("NON_LINEAR");
    parser.add_enum_option("internal_merging",
                           internal_merging,
                           "choose the order in which to merge the set of "
                           "abstractions to be merged (only useful with "
                           "MERGE_FOR_ATOMIC): "
                           "linear (obvious), "
                           "non linear, which means to first merge every cycle, "
                           "and then the resulting intermediate abstractions.",
                           "LINEAR");

    // Options for GraphCreator
    parser.add_option<bool>("build_stabilized_pdg", "build an abstraction "
                            "stabilized pdb, which results in bliss searching "
                            "for local symmetries only", "False");
    parser.add_option<bool>("debug_graph_creator", "produce dot readable output "
                            "from the graph generating methods", "false");

    Options options = parser.parse();
    if (options.get<int>("bliss_call_time_limit")
            && options.get<int>("bliss_total_time_budget")) {
        cerr << "Please only specify bliss_call_time_limite or "
                "bliss_total_time_budget but not both" << endl;
        exit_with(EXIT_INPUT_ERROR);
    }
    if (options.get_enum("symmetries_for_shrinking") == 0
            && options.get_enum("symmetries_for_merging") == 0) {
        cerr << "Please use symmetries at least for shrinking or merging, "
                "otherwise use merge_dfp instead." << endl;
        exit_with(EXIT_INPUT_ERROR);
    }
    if (parser.dry_run())
        return 0;
    else
        return new MergeSymmetriesLinear(options);
}

static Plugin<MergeStrategy> _plugin("merge_symmetries_linear", _parse);
