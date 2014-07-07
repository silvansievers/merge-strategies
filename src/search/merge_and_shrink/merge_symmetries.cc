#include "merge_symmetries.h"

#include "symmetries/symmetries.h"

#include "../plugin.h"

#include <limits>

using namespace std;

MergeSymmetries::MergeSymmetries(const Options &options_)
    : MergeDFP(),
      options(options_),
      max_symmetry_iterations(options.get<int>("max_symmetry_iterations")),
      internal_merging(InternalMerging(options.get_enum("internal_merging"))),
      started_merging_for_symmetries(false),
      number_of_applied_symmetries(0),
      atomic_symmetries(0),
      binary_symmetries(0),
      other_symmetries(0),
      iteration_counter(0) {
}

void MergeSymmetries::dump_statistics() const {
    if (iteration_counter == 1) {
        cout << "First iteration: atomic symmetries: " << atomic_symmetries << endl;
        cout << "First iteration: binary symmetries: " << binary_symmetries << endl;
        cout << "First iteration: other symmetries: " << other_symmetries << endl;
    } else if (remaining_merges == 1) {
        cout << "Total: atomic symmetries: " << atomic_symmetries << endl;
        cout << "Total: binary symmetries: " << binary_symmetries << endl;
        cout << "Total: other symmetries: " << other_symmetries << endl;
    }
}

pair<int, int> MergeSymmetries::get_next(vector<Abstraction *> &all_abstractions) {
    assert(!done());
    ++iteration_counter;

    if (iteration_counter <= max_symmetry_iterations && abs_to_merge.empty()) {
        Symmetries symmetries(options);
        if (started_merging_for_symmetries) {
            // TODO: can we somehow make sure that if we were merging in order
            // to apply a symmetry and no shrinking happened, then we indeed
            // applied an atomic symmetry in the line above?
            started_merging_for_symmetries = false;
        }
        pair<int, int> stats = symmetries.find_and_apply_symmetries(all_abstractions, abs_to_merge);
        number_of_applied_symmetries += stats.first;
        remaining_merges -= stats.second;
        cout << "Number of applied symmetries: " << number_of_applied_symmetries << endl;
        atomic_symmetries += symmetries.get_atomic_symmetries();
        binary_symmetries += symmetries.get_binary_symmetries();
        other_symmetries += symmetries.get_other_symmetries();
        if (abs_to_merge.size() == 1) {
            cerr << "Atomic symmetry, not applied!" << endl;
            exit_with(EXIT_CRITICAL_ERROR);
        } else if (abs_to_merge.size() == 0) {
            cout << "No symmetries for merging found." << endl;
        } else {
            cout << "Merging next: ";
            for (vector<int>::iterator it = abs_to_merge.begin(); it != abs_to_merge.end(); ++it) {
                cout << *it << " ";
            }
            cout << endl;
        }
    }

    dump_statistics();

    if (abs_to_merge.empty() || iteration_counter > max_symmetry_iterations) {
        return MergeDFP::get_next(all_abstractions);
    }

    int first;
    if (!started_merging_for_symmetries) {
        assert(abs_to_merge.size() > 1);
        started_merging_for_symmetries = true;
        first = *abs_to_merge.begin();
        abs_to_merge.erase(abs_to_merge.begin());
    } else {
        first = all_abstractions.size() - 1;
    }
    assert(!abs_to_merge.empty());
    int second = *abs_to_merge.begin();
    abs_to_merge.erase(abs_to_merge.begin());

    --remaining_merges;
    return make_pair(first, second);
}

string MergeSymmetries::name() const {
    return "symmetries";
}

static MergeStrategy *_parse(OptionParser &parser) {
    parser.add_option<bool>("debug_graph_creator", "produce dot readable output "
                            "from the graph generating methods", "false");
    parser.add_option<int>("max_symmetry_iterations", "number of iteration up "
                           "to which symmetries should be searched for and "
                           "applied.", "infinity");
    vector<string> symmetries_for_shrinking;
    symmetries_for_shrinking.push_back("ATOMIC");
    symmetries_for_shrinking.push_back("LOCAL");
    symmetries_for_shrinking.push_back("NONE");
    parser.add_enum_option("symmetries_for_shrinking",
                           symmetries_for_shrinking,
                           "choose the type of symmetries used for shrinking: "
                           "only atomic symmetries, "
                           "local symmetries, "
                           "only use for merging, no shrinking.",
                           "ATOMIC");
    vector<string> internal_merging;
    internal_merging.push_back("LINEAR");
    internal_merging.push_back("NON_LINEAR");
    parser.add_enum_option("internal_merging",
                           internal_merging,
                           "choose how the set of abstractions that must be "
                           "merged for symmetries is merged: "
                           "linearly, resulting in one large abstractions, "
                           "non linearly, resulting in one composite abstraction "
                           "for every previous cycle of the chosen symmetry.",
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
