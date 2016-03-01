#include "merge_random_linear.h"

#include "factored_transition_system.h"

#include "../option_parser.h"
#include "../plugin.h"
#include "../task_proxy.h"

#include "../utils/collections.h"
#include "../utils/memory.h"
#include "../utils/rng.h"

#include <cassert>
#include <iostream>

using namespace std;

namespace merge_and_shrink {
MergeRandomLinear::MergeRandomLinear(const Options &options)
    : MergeStrategy(),
      random_seed(options.get<int>("random_seed")),
      need_first_index(true) {
}

void MergeRandomLinear::initialize(const shared_ptr<AbstractTask> task) {
    MergeStrategy::initialize(task);
    TaskProxy task_proxy(*task);
    int num_variables = task_proxy.get_variables().size();
    vector<int> variables(num_variables, -1);
    iota(variables.begin(), variables.end(), 0);

    utils::RandomNumberGenerator rng(random_seed);
    randomized_variable_order.reserve(num_variables);
    for (int i = 0; i < num_variables; ++i) {
        int random_index = rng(variables.size());
        randomized_variable_order.push_back(variables[random_index]);
        variables.erase(variables.begin() + random_index);
    }
}

pair<int, int> MergeRandomLinear::get_next(
    FactoredTransitionSystem &fts) {
    assert(initialized());
    assert(!done());

    int next_index1;
    if (need_first_index) {
        need_first_index = false;
        next_index1 = randomized_variable_order.front();
        randomized_variable_order.erase(randomized_variable_order.begin());
        cout << "First variable: " << next_index1 << endl;
    } else {
        // The most recent composite transition system is appended at the end
        int num_transition_systems = fts.get_size();
        next_index1 = num_transition_systems - 1;
    }
    int next_index2 = randomized_variable_order.front();
    randomized_variable_order.erase(randomized_variable_order.begin());
    cout << "Next variable: " << next_index2 << endl;
    assert(fts.is_active(next_index1));
    assert(fts.is_active(next_index2));
    --remaining_merges;
    return make_pair(next_index1, next_index2);
}

void MergeRandomLinear::dump_strategy_specific_options() const {
    cout << "random seed: " << random_seed << endl;
}

string MergeRandomLinear::name() const {
    return "random linear";
}

static shared_ptr<MergeStrategy>_parse(OptionParser &parser) {
    parser.document_synopsis(
        "Random linear merge strategy.",
        "This merge strategy randomly computes a variable order for merging.");
    parser.add_option<int>("random_seed", "random seed", "2015");

    Options opts = parser.parse();
    if (parser.dry_run())
        return nullptr;
    else
        return make_shared<MergeRandomLinear>(opts);
}

static PluginShared<MergeStrategy> _plugin("merge_random_linear", _parse);
}
