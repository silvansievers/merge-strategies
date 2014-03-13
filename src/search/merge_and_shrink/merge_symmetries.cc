#include "merge_symmetries.h"

#include "symmetries/symmetries.h"

#include "../option_parser.h"
#include "../plugin.h"

using namespace std;

MergeSymmetries::MergeSymmetries(const Options &options_)
    : MergeDFP(), options(options_), symmetries(0), index_of_composite_abs(-1) {
}

bool MergeSymmetries::done() const {
    return MergeDFP::done();
}

pair<int, int> MergeSymmetries::get_next(const vector<Abstraction *> &all_abstractions) {
    // TODO: how and when do we invalidate the symmetries object
    if (abs_to_merge.empty()) {
        bool found_symmetry = symmetries->find_and_apply_atomar_symmetries(all_abstractions);
        if (found_symmetry) {
            cout << "Found and applied atomar symmetries." << endl;
        } else {
            cout << "No atomar symmetries found." << endl;
        }
        found_symmetry = symmetries->find_to_be_merged_abstractions(all_abstractions, abs_to_merge);
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
    int second;
    if (index_of_composite_abs == -1) {
        first = *abs_to_merge.begin();
        abs_to_merge.erase(abs_to_merge.begin());
    } else {
        first = index_of_composite_abs;
    }
    second = *abs_to_merge.begin();
    abs_to_merge.erase(abs_to_merge.begin());
    return make_pair(first, second);
}

string MergeSymmetries::name() const {
    return "symmetries";
}

void MergeSymmetries::initialize(const Labels *labels) {
    symmetries = new Symmetries(options, labels);
}

static MergeStrategy *_parse(OptionParser &parser) {
    parser.add_option<bool>("debug_graph_creator", "produce dot readable output "
                            "from the graph generating methods", "false");
    parser.add_option<int>("version", "debug application of atomar symmetries", "1");
    Options options = parser.parse();
    if (!parser.dry_run())
        return new MergeSymmetries(options);
    else
        return 0;
}

static Plugin<MergeStrategy> _plugin("merge_symmetries", _parse);
