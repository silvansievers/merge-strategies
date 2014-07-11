#include "merge_symmetries.h"

#include "symmetries/symmetries.h"

#include "../plugin.h"

#include <limits>

using namespace std;

MergeSymmetries::MergeSymmetries(const Options &options_)
    : MergeDFP(),
      options(options_),
      max_iterations(options.get<int>("max_iterations")),
      number_of_applied_symmetries(0),
      iteration_counter(0) {
}

pair<int, int> MergeSymmetries::get_next(vector<Abstraction *> &all_abstractions) {
    assert(!done());
    ++iteration_counter;

    if (iteration_counter <= max_iterations && merge_order.empty()) {
        Symmetries symmetries(options);
        bool applied_symmetries =
                symmetries.find_and_apply_symmetries(all_abstractions, merge_order);
        if (applied_symmetries)
            ++number_of_applied_symmetries;
        cout << "Number of applied symmetries: " << number_of_applied_symmetries << endl;
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
    parser.add_option<int>("max_iterations", "number of merge-and-shrink "
                           "iterations up to which symmetries should be computed."
                           "infinity");
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
