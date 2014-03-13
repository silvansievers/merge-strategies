#include "symmetries_merge_strategy.h"

SymmetriesMergeStrategy::SymmetriesMergeStrategy(const Options &opts)
    : NonLinearMergeStrategy(opts) {
}

SymmetriesMergeStrategy::done() const {
    return NonLinearMergeStrategy::done();
}

void NonLinearMergeStrategy::get_next(const std::vector<Abstraction *> &all_abstractions, pair<int, int> &next_indices) {
}

void NonLinearMergeStrategy::dump_strategy_specific_options() const {
    cout << "Non linear merge strategy type: ";
    switch (non_linear_merge_strategy_type) {
    case DFP:
        cout << "DFP";
        break;
    default:
        ABORT("Unknown merge strategy.");
    }
    cout << endl;
}

string NonLinearMergeStrategy::name() const {
    return "non linear";
}

static MergeStrategy *_parse(OptionParser &parser) {
    vector<string> merge_strategies;
    //TODO: it's a bit dangerous that the merge strategies here
    // have to be specified exactly in the same order
    // as in the enum definition. Try to find a way around this,
    // or at least raise an error when the order is wrong.
    merge_strategies.push_back("DFP");
    parser.add_enum_option("type", merge_strategies,
                           "non linear merge strategy",
                           "DFP");

    Options opts = parser.parse();
    if (parser.help_mode())
        return 0;
    if (!parser.dry_run())
        return new NonLinearMergeStrategy(opts);
    else
        return 0;
}

static Plugin<MergeStrategy> _plugin("merge_non_linear", _parse);
