#include "merge_symmetries.h"

#include "symmetries/symmetries.h"

#include "../plugin.h"

#include <limits>

using namespace std;

MergeSymmetries::MergeSymmetries(const Options &options_)
    : MergeLinear(options_),
      options(options_),
      started_merging_for_symmetries(false),
      number_of_applied_symmetries(0),
      atomic_symmetries(0),
      binary_symmetries(0),
      other_symmetries(0),
      iteration_counter(0),
      max_symmetry_iterations(options.get<int>("max_symmetry_iterations")) {
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

bool MergeSymmetries::done() const {
    return MergeLinear::done();
}

pair<int, int> MergeSymmetries::get_next(const vector<Abstraction *> &all_abstractions) {
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
        number_of_applied_symmetries += symmetries.find_and_apply_symmetries(all_abstractions, abs_to_merge);
        cout << "Number of applied symmetries: " << number_of_applied_symmetries << endl;
        atomic_symmetries += symmetries.get_atomic_symmetries();
        binary_symmetries += symmetries.get_binary_symmetries();
        other_symmetries += symmetries.get_other_symmetries();
        if (abs_to_merge.size() == 1) {
            cerr << "Atomic symmetry, not applied!" << endl;
            exit_with(EXIT_CRITICAL_ERROR);
        } else if (abs_to_merge.size() == 0) {
            cout << "No symmetries for merging found." << endl;
        }
    }

    dump_statistics();

    if (abs_to_merge.empty() || iteration_counter > max_symmetry_iterations) {
        return MergeLinear::get_next(all_abstractions);
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
    return "linear";
}

static MergeStrategy *_parse(OptionParser &parser) {
    parser.add_option<bool>("debug_graph_creator", "produce dot readable output "
                            "from the graph generating methods", "false");
    parser.add_option<int>("max_symmetry_iterations", "number of iteration up "
                           "to which symmetries should be searched for and "
                           "applied.", "infinity");
    vector<string> type_of_symmetries;
    type_of_symmetries.push_back("ATOMIC");
    type_of_symmetries.push_back("LOCAL");
    type_of_symmetries.push_back("ANY");
    parser.add_enum_option("type_of_symmetries", type_of_symmetries,
                           "typo of symmetrie: only atomic symmetries, "
                           "local and atomic symmetries, "
                           "any kind of symmetry.", "ATOMIC");
    parser.add_option<bool>("build_stabilized_pdg", "build an abstraction "
                            "stabilized pdb, which results in bliss searching "
                            "for local symmetries only", "False");
    vector<string> merge_strategies;
    //TODO: it's a bit dangerous that the merge strategies here
    // have to be specified exactly in the same order
    // as in the enum definition. Try to find a way around this,
    // or at least raise an error when the order is wrong.
    merge_strategies.push_back("CG_GOAL_LEVEL");
    merge_strategies.push_back("CG_GOAL_RANDOM");
    merge_strategies.push_back("GOAL_CG_LEVEL");
    merge_strategies.push_back("RANDOM");
    merge_strategies.push_back("LEVEL");
    merge_strategies.push_back("REVERSE_LEVEL");
    parser.add_enum_option("variable_order", merge_strategies,
                           "the order in which atomic abstractions are merged",
                           "CG_GOAL_LEVEL");

    Options options = parser.parse();
    if (parser.dry_run())
        return 0;
    else
        return new MergeSymmetries(options);
}

static Plugin<MergeStrategy> _plugin("merge_symmetries", _parse);
