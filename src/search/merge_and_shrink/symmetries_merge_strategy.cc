#include "symmetries_merge_strategy.h"

#include <iostream>

using namespace std;

SymmetriesMergeStrategy::SymmetriesMergeStrategy(const Options &opts)
    : NonLinearMergeStrategy(opts) {
}

bool SymmetriesMergeStrategy::done() const {
    return NonLinearMergeStrategy::done();
}

pair<int, int> SymmetriesMergeStrategy::get_next(const vector<Abstraction *> &all_abstractions) {
    return make_pair(all_abstractions.size(), 0);
}

void SymmetriesMergeStrategy::dump_strategy_specific_options() const {
    cout << "Symmetries merge strategy with underlying merge strategy:" << endl;
    SymmetriesMergeStrategy::dump_strategy_specific_options();
}

string SymmetriesMergeStrategy::name() const {
    return "symmetries";
}

//static MergeStrategy *_parse(OptionParser &parser) {
//    vector<string> merge_strategies;
//    //TODO: it's a bit dangerous that the merge strategies here
//    // have to be specified exactly in the same order
//    // as in the enum definition. Try to find a way around this,
//    // or at least raise an error when the order is wrong.
//    merge_strategies.push_back("DFP");
//    parser.add_enum_option("type", merge_strategies,
//                           "non linear merge strategy",
//                           "DFP");

//    Options opts = parser.parse();
//    if (parser.help_mode())
//        return 0;
//    if (!parser.dry_run())
//        return new NonLinearMergeStrategy(opts);
//    else
//        return 0;
//}

//static Plugin<MergeStrategy> _plugin("merge_non_linear", _parse);
