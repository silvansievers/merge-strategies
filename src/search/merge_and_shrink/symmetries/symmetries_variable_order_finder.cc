#include "symmetries_variable_order_finder.h"

#include "../abstraction.h"
#include "../../utilities.h"

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <vector>

using namespace std;

SymmetriesVariableOrderFinder::SymmetriesVariableOrderFinder(MergeStrategy merge_strategy,
                                                             bool is_first,
                                                             bool debug_graph_creator,
                                                             int version,
                                                             const vector<Abstraction *> &atomic_abstractions)
    : VariableOrderFinder(merge_strategy, is_first),
      symmetries(debug_graph_creator, version) {
    found_symmetry = symmetries.find_and_apply_atomar_symmetries(atomic_abstractions);
    if (found_symmetry) {
        cout << "Found and applied atomar symmetries on atomic abstractions." << endl;
    } else {
        cout << "No atomar symmetries found for atomic abstractions." << endl;
    }
    found_symmetry = false;
    // initialize by searching for symmetries.
    compute_symmetries(atomic_abstractions);
}

void SymmetriesVariableOrderFinder::set_used_var(int var_no) {
    for (size_t i = 0; i < remaining_vars.size(); ++i) {
        if (remaining_vars[i] == var_no) {
            select_next(i, var_no);
            break;
        }
        //if (i == remaining_vars.size() - 1)
        //    assert(false);
        // NOTE: it may happen that this function is called although the variable
        // has already been marked as "used" in the variable order finder, because
        // we may start merging for symmetries from the same index several times
    }
}

void SymmetriesVariableOrderFinder::compute_symmetries(const vector<Abstraction *> &abstractions) {
    assert(!found_symmetry);
    assert(abs_to_merge.empty());
    found_symmetry = symmetries.find_to_be_merged_abstractions(abstractions, abs_to_merge);
    if (found_symmetry) {
        if (abs_to_merge.size() <= 1) {
            cout << "Found an atomar symmetry which will be immediately applied." << endl;
        } else {
            cout << "Found symmetries. Will merge accordingly in future iteration(s)." << endl;
        }
    } else {
        cout << "No symmetries found, use chosen standard M&S options." << endl;
    }
}

int SymmetriesVariableOrderFinder::get_next_to_merge() {
    assert(found_symmetry);
    assert(!abs_to_merge.empty());
    int var_no = -1;
    cout << "Need to merge more abstractions before application of symmetries." << endl;
    var_no = *abs_to_merge.begin();
    abs_to_merge.erase(abs_to_merge.begin());
    set_used_var(var_no);
    return var_no;
}

void SymmetriesVariableOrderFinder::apply_symmetries(const vector<Abstraction *> &abstractions) {
    assert(found_symmetry);
    assert(abs_to_merge.empty());
    cout << "Apply atomar symmetry which now has to be applicable." << endl;
    found_symmetry = symmetries.find_and_apply_atomar_symmetries(abstractions);
    assert(found_symmetry);
    found_symmetry = false;
}
