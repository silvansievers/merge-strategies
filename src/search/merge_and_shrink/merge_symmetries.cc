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
    if (abs_to_merge.empty()) {
        Symmetries symmetries(options);
        /* Some thoughts: can we combine find and apply symmetries and finde to be merged abstractions
         * into one method? e.g., here we have something like
         * while (find to be merged abstractions) {
         *     if !abs_to_merge.empty()
         *         break // then we need to actually merge.
         * }
         */
        bool found_symmetry = symmetries.find_and_apply_atomar_symmetries(all_abstractions);
        if (started_merging_for_symmetries) {
            // we were merging abstractions to allow symmetries application
            // TODO: what if shrinking happens?
            assert(found_symmetry);
            started_merging_for_symmetries = false;
        }
        if (found_symmetry) {
            cout << "Found and applied atomar symmetries." << endl;
        } else {
            cout << "No atomar symmetries found." << endl;
        }
        found_symmetry = symmetries.find_to_be_merged_abstractions(all_abstractions, abs_to_merge);
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
    return make_pair(first, second);
}

string MergeSymmetries::name() const {
    return "symmetries";
}

static MergeStrategy *_parse(OptionParser &parser) {
    parser.add_option<bool>("debug_graph_creator", "produce dot readable output "
                            "from the graph generating methods", "false");
    parser.add_option<int>("version", "debug application of atomar symmetries", "1");

    Options options = parser.parse();
    if (parser.dry_run())
        return 0;
    else
        return new MergeSymmetries(options);
}

static Plugin<MergeStrategy> _plugin("merge_symmetries", _parse);
