#include "merge_predefined.h"

#include "../option_parser.h"
#include "../plugin.h"
#include "../utilities.h"

#include <cassert>
#include <iostream>

using namespace std;

MergePredefined::MergePredefined(const Options &options)
    : MergeStrategy(),
      merge_order(options.get_list<vector<int>>("merge_order")) {
}

void MergePredefined::initialize(const shared_ptr<AbstractTask> task) {
    MergeStrategy::initialize(task);
    if (static_cast<int>(merge_order.size()) != remaining_merges) {
        cout << remaining_merges << endl;
        cerr << "Invalid size of merge order" << endl;
        exit_with(EXIT_INPUT_ERROR);
    }
}

pair<int, int> MergePredefined::get_next(
    const vector<TransitionSystem *> &all_transition_systems) {
    assert(!merge_order.empty());
    const vector<int> &next_pair = merge_order.front();
    assert(next_pair.size() == 2);
    int next_index1 = next_pair[0];
    int next_index2 = next_pair[1];
    assert(all_transition_systems[next_index1]);
    assert(all_transition_systems[next_index2]);
    merge_order.erase(merge_order.begin());

    --remaining_merges;
    cout << "Next pair of indices: (" << next_index1 << ", "
         << next_index2 << ")" << endl;
    return make_pair(next_index1, next_index2);
}

void MergePredefined::dump_strategy_specific_options() const {
    cout << "merge order: " << merge_order << endl;
}

string MergePredefined::name() const {
    return "predefined";
}

static shared_ptr<MergeStrategy>_parse(OptionParser &parser) {
    parser.document_synopsis(
        "Random merge strategy.",
        "This merge strategy randomly selects the two next transition systems"
        "to merge.");
    parser.add_list_option<vector<int>>(
        "merge_order",
        "predefined merge order",
        "");

    Options options = parser.parse();
    if (parser.dry_run()) {
        vector<vector<int>> merge_order = options.get_list<vector<int>>("merge_order");
        if (merge_order.empty()) {
            cerr << "Got empty merge order, aborting" << endl;
            exit_with(EXIT_INPUT_ERROR);
        }
        for (const vector<int> &pair : merge_order) {
            if (pair.size() != 2) {
                cerr << "Every element in the list merge_order must contain "
                        "exactly two elements!" << endl;
                cout << pair << endl;
                exit_with(EXIT_INPUT_ERROR);
            }
        }
        return nullptr;
    } else {
        return make_shared<MergePredefined>(options);
    }
}

static PluginShared<MergeStrategy> _plugin("merge_predefined", _parse);
