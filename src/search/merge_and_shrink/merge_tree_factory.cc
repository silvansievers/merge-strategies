#include "merge_tree_factory.h"

#include "merge_tree.h"

#include "../options/option_parser.h"
#include "../options/plugin.h"

#include "../utils/rng_options.h"
#include "../utils/system.h"

#include <iostream>

using namespace std;

namespace merge_and_shrink {
MergeTreeFactory::MergeTreeFactory(const options::Options &options)
    : update_option(static_cast<UpdateOption>(options.get_enum("update_option"))) {
    rng = utils::parse_rng_from_options(options);
}

void MergeTreeFactory::dump_options() const {
    cout << "Merge tree options: " << endl;
    cout << "Type: " << name() << endl;
    cout << "Update option: ";
    switch (update_option) {
    case UpdateOption::USE_FIRST:
        cout << "use first";
        break;
    case UpdateOption::USE_SECOND:
        cout << "use second";
        break;
    case UpdateOption::USE_RANDOM:
        cout << "use random";
        break;
    case UpdateOption::USE_GREATER_INDEX:
        cout << "use greater index";
        break;
    case UpdateOption::USE_SMALLER_INDEX:
        cout << "use smaller index";
        break;
    }
    dump_tree_specific_options();
}

unique_ptr<MergeTree> MergeTreeFactory::compute_merge_tree(
    shared_ptr<AbstractTask>,
    FactoredTransitionSystem &,
    const std::vector<int> &) {
    cerr << "This merge tree does not support being computed on a subset "
            "of indices for a given factored transition system!" << endl;
    utils::exit_with(utils::ExitCode::CRITICAL_ERROR);
}

void MergeTreeFactory::add_options_to_parser(options::OptionParser &parser) {
    utils::add_rng_options(parser);
    vector<string> update_option;
    update_option.push_back("use_first");
    update_option.push_back("use_second");
    update_option.push_back("use_random");
    update_option.push_back("use_greater_index");
    update_option.push_back("use_smaller_index");
    parser.add_enum_option(
        "update_option",
        update_option,
        "When the merge tree is used within another merge strategy, how"
        "should it be updated when a merge different to a merge from the "
        "tree is performed: choose among use_first, use_second, "
        "use_random, use_greater_index, and use_smaller_index to let the left "
        "(right, random, greater, smaller) leaf node in the tree survive and "
        "remove the other leaf node without replacement. ",
        "use_random");

}

static options::PluginTypePlugin<MergeTreeFactory> _type_plugin(
    "Merge Tree",
    "This page describes the available merge trees that can be used to "
    "precompute a merge strategy, either for the entire task or a given "
    "subset of transition systems of a given factored transition system.\n"
    "Merge trees are used in the merge strategies of type 'precomputed'.");
}
