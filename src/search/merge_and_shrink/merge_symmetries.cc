#include "merge_symmetries.h"

#include "../option_parser.h"
#include "../plugin.h"

using namespace std;

MergeSymmetries::MergeSymmetries()
    : MergeDFP() {
}

bool MergeSymmetries::done() const {
    return MergeDFP::done();
}

pair<int, int> MergeSymmetries::get_next(const Labels *,
                                         const vector<Abstraction *> &all_abstractions) {
    return make_pair(all_abstractions.size(), 0);
}

string MergeSymmetries::name() const {
    return "symmetries";
}

static MergeStrategy *_parse(OptionParser &parser) {
    if (!parser.dry_run())
        return new MergeSymmetries();
    else
        return 0;
}

static Plugin<MergeStrategy> _plugin("merge_symmetries", _parse);
