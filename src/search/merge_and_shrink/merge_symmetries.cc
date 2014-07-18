#include "merge_symmetries.h"

#include "symmetries/symmetries.h"

#include "../plugin.h"

#include <algorithm>
#include <iomanip>

using namespace std;

MergeSymmetries::MergeSymmetries(const Options &options_)
    : MergeDFP(),
      options(options_),
      number_of_applied_symmetries(0),
      bliss_limit_reached(false) {
}

void MergeSymmetries::dump_statistics() {
    if (!bliss_times.empty()) {
        sort(bliss_times.begin(), bliss_times.end());
        double summed_up_bliss_times = 0;
        size_t total_bliss_calls = bliss_times.size();
        cout << "total bliss calls: " << total_bliss_calls << endl;
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

pair<int, int> MergeSymmetries::get_next(vector<Abstraction *> &all_abstractions) {
    assert(!done());

    if (!bliss_limit_reached && merge_order.empty()) {
        Symmetries symmetries(options);
        bool applied_symmetries =
                symmetries.find_and_apply_symmetries(all_abstractions, merge_order);
        if (applied_symmetries)
            ++number_of_applied_symmetries;
        if (symmetries.is_bliss_limit_reached()) {
            bliss_limit_reached = true;
        }
        bliss_times.push_back(symmetries.get_bliss_time());
        cout << "Number of applied symmetries: " << number_of_applied_symmetries << endl;
    }

    if (remaining_merges == 1) {
        dump_statistics();
    }

    if (merge_order.empty()) {
        return MergeDFP::get_next(all_abstractions);
    }

    pair<int, int> next_merge = merge_order.front();
    if (!all_abstractions[next_merge.first] || !all_abstractions[next_merge.second]) {
        cerr << "Problem with the merge strategy: invalid indices" << endl;
        exit_with(EXIT_CRITICAL_ERROR);
    }
    merge_order.erase(merge_order.begin());

    --remaining_merges;
    return next_merge;
}

string MergeSymmetries::name() const {
    return "symmetries";
}

static MergeStrategy *_parse(OptionParser &parser) {
    parser.add_option<bool>("debug_graph_creator", "produce dot readable output "
                            "from the graph generating methods", "false");
    parser.add_option<int>("bliss_time_limit", "time in seconds one bliss "
                           "run is allowed to last at most (0 means no limit)",
                           "0");
    vector<string> symmetries_for_shrinking;
    symmetries_for_shrinking.push_back("NO_SHRINKING");
    symmetries_for_shrinking.push_back("ATOMIC");
    symmetries_for_shrinking.push_back("LOCAL");
    symmetries_for_shrinking.push_back("LOCAL_ONE_ABS");
    symmetries_for_shrinking.push_back("LOCAL_ONE_SYM");
    symmetries_for_shrinking.push_back("LOCAL_ONE_ABS_SYM");
    parser.add_enum_option("symmetries_for_shrinking",
                           symmetries_for_shrinking,
                           "choose the type of symmetries used for shrinking: "
                           "no shrinking, "
                           "only atomic symmetries, "
                           "local symmetries.",
                           "NO_SHRINKING");
    vector<string> symmetries_for_merging;
    symmetries_for_merging.push_back("NO_MERGING");
    symmetries_for_merging.push_back("LEAST_OVERALL_AFFECTED");
    symmetries_for_merging.push_back("MOST_OVERALL_AFFECTED");
    symmetries_for_merging.push_back("LEAST_MAPPED");
    symmetries_for_merging.push_back("MOST_MAPPED");
    parser.add_enum_option("symmetries_for_merging",
                           symmetries_for_merging,
                           "choose the type of symmetries used for merging: "
                           "the smallest or the largest in number of abstractions "
                           "that are affected or mapped.",
                           "LEAST_OVERALL_AFFECTED");
    vector<string> internal_merging;
    internal_merging.push_back("LINEAR");
    internal_merging.push_back("NON_LINEAR");
    internal_merging.push_back("NON_LINEAR_INCOMPLETE");
    parser.add_enum_option("internal_merging",
                           internal_merging,
                           "choose how the set of abstractions that must be "
                           "merged for symmetries is merged: "
                           "linearly, over the entire set of abstractions, "
                           "non linearly, which means merging every cycle, "
                           "then linearly all remaining and resultin "
                           "abstractions.",
                           "LINEAR");
    parser.add_option<bool>("build_stabilized_pdg", "build an abstraction "
                            "stabilized pdb, which results in bliss searching "
                            "for local symmetries only", "False");

    Options options = parser.parse();
    if (parser.dry_run())
        return 0;
    else
        return new MergeSymmetries(options);
}

static Plugin<MergeStrategy> _plugin("merge_symmetries", _parse);
