#include "merge_symmetries.h"

#include "symmetries/symmetries.h"

#include "../plugin.h"

using namespace std;

MergeSymmetries::MergeSymmetries(const Options &options_)
    : MergeDFP(), options(options_), started_merging_for_symmetries(false) {
}

bool MergeSymmetries::done() const {
    return MergeDFP::done();
}

pair<int, int> MergeSymmetries::get_next(const vector<Abstraction *> &all_abstractions) {
    assert(!done());

    if (abs_to_merge.empty()) {
        Symmetries symmetries(options);
        /* Some thoughts: can we combine find and apply symmetries and find
         * to be merged abstractions into one method? e.g., we could have
         * something like at this place:
         * while (find to be merged abstractions) {
         *     if !abs_to_merge.empty()
         *         break // then we need to actually merge.
         * }
         */
        // We must assert that all abstractions distances have been computed
        // because of the nasty possible side effect of pruning irrelevant
        // states. Alternativeley, we could compute distances here, *before*
        // searching for symmetries.
        for (size_t i = 0; i < all_abstractions.size(); ++i) {
            if (all_abstractions[i])
                assert(all_abstractions[i]->are_distances_computed());
        }
        //symmetries.find_and_apply_atomic_symmetries(all_abstractions);
        if (started_merging_for_symmetries) {
            // TODO: can we somehow make sure that if we were merging in order
            // to apply a symmetry and no shrinking happened, then we indeed
            // applied an atomic symmetry in the line above?
            started_merging_for_symmetries = false;
        }
        bool found_symmetry = symmetries.find_to_be_merged_abstractions(all_abstractions, abs_to_merge);
        if (found_symmetry) {
            assert(abs_to_merge.size() > 1);
        } else {
            cout << "No symmetries found at all." << endl;
        }
    }

    if (abs_to_merge.empty()) {
        return MergeDFP::get_next(all_abstractions);
    }

    int first;
    if (!started_merging_for_symmetries) {
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
    parser.add_option<int>("version", "apply all atomic symmetries at the same time (1) or separatedly (0)", "1");

    Options options = parser.parse();
    if (parser.dry_run())
        return 0;
    else
        return new MergeSymmetries(options);
}

static Plugin<MergeStrategy> _plugin("merge_symmetries", _parse);
