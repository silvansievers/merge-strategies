#include "symmetries.h"

#include "permutation.h"
#include "scc.h"

#include "../abstraction.h"

#include "../../globals.h"
#include "../../option_parser.h"
#include "../../timer.h"
#include "../../utilities.h"

#include <cassert>
#include <iostream>
#include <limits>

// TODO: copied and renamed from shrink_strategy.h
typedef __gnu_cxx::slist<AbstractStateRef> EquivClass;
typedef std::vector<EquivClass> EquivRel;

Symmetries::Symmetries(const Options &options)
    : gc(options.get<bool>("debug_graph_creator")),
      version(options.get<int>("version")) {
}

bool Symmetries::is_atomic_generator(const vector<Abstraction *> abstractions, int gen_index) const {
    unsigned int num_abstractions = get_permutations_wrapper().num_abstractions;
    // Check if an abstraction is not mapped to itself
    for (unsigned int index = 0; index < num_abstractions; ++index) {
        if (abstractions[index] == 0)
            continue;
        unsigned int to_index = get_generator(gen_index)->get_value(index);
        if (index == to_index)
            continue;
        return false;
    }
    // Check if more than one abstraction is affected non-trivially
    set<int> affected_abs;
    for (unsigned int index = num_abstractions; index < get_permutations_wrapper().num_abs_and_states; ++index) {
        int abs_index = get_permutations_wrapper().get_var_by_index(index);
        if (abstractions[abs_index] == 0)
            continue;
        unsigned int to_index = get_generator(gen_index)->get_value(index);
        if (index == to_index)
            continue;
        if (affected_abs.insert(abs_index).second) {
            cout << "abstraction " << abstractions[abs_index]->description() << " non-trivially affected" << endl;
            if (affected_abs.size() > 1) {
                return false;
            }
        }
    }
    return true;
}

void Symmetries::find_and_apply_atomic_symmetries(const vector<Abstraction *> &abstractions) {
    cout << "Looking for atomic symmetries:" << endl;
    while (true) {
        vector<set<int> > non_trivially_affected_abstractions;
        vector<bool> atomic_symmetries;
        vector<vector<int> > atomic_symmetries_by_affected_abs;
        bool found_symmetry = find_symmetries(abstractions, non_trivially_affected_abstractions, atomic_symmetries, atomic_symmetries_by_affected_abs, true);
        if (!found_symmetry) {
            return;
        }
        if (version == 1 || version == 0) {
            int most_affected_abs_index = -1;
            int most_affected_abs_size = 0;
            for (size_t abs_index = 0; abs_index < abstractions.size(); ++abs_index) {
                if (abstractions[abs_index]) {
                    int no_atomic_symmetries_affected = atomic_symmetries_by_affected_abs[abs_index].size();
                    if (no_atomic_symmetries_affected > most_affected_abs_size) {
                        most_affected_abs_size = no_atomic_symmetries_affected;
                        most_affected_abs_index = abs_index;
                    }
                }
            }
            assert(most_affected_abs_index != -1);
            assert(most_affected_abs_size > 0);
            // We apply all symmetries which affect one abstraction.
            apply_symmetries(abstractions, atomic_symmetries_by_affected_abs[most_affected_abs_index]);
        }
    }
}

bool Symmetries::find_atomic_symmetries(const vector<Abstraction *>& abstractions,
                                        vector<vector<int> > &atomic_symmetries_by_affected_abs) {
    cout << "Computing generators for atomic symmetries" << endl;
    gc.compute_generators(abstractions, true);
    if (get_num_generators() == 0) {
        cout << "No generators found! Done searching for atomic symmetries. [t=" << g_timer << "]" << endl;
        return false;
    }

    int num_abstractions = get_permutations_wrapper().num_abstractions;
    atomic_symmetries_by_affected_abs.resize(abstractions.size(), vector<int>());
    for (int gen_index = 0; gen_index < get_num_generators(); ++gen_index) {
        for (unsigned int index = 0; index < num_abstractions; ++index) {
            if (abstractions[index] == 0) {
                unsigned int to_index = get_generator(gen_index)->get_value(index);
                if (index != to_index) {
                    cerr << "non atomic symmetry found" << endl;
                    exit_with(EXIT_CRITICAL_ERROR);
                }
            }
        }

        set<int> affected_abstractions;
        for (unsigned int index = num_abstractions; index < get_permutations_wrapper().num_abs_and_states; ++index) {
            int abs_index = get_permutations_wrapper().get_var_by_index(index);
            // TODO: can we assert(abstractions[abs_index]), because empty
            // abstractions have no states?
            if (abstractions[abs_index]) {
                unsigned int to_index = get_generator(gen_index)->get_value(index);
                if (index != to_index) {
                    // TODO: check if generator already affects another abstraction
                    if (affected_abstractions.insert(abs_index).second) {
                        cout << "abstraction " << abstractions[abs_index]->description() << " non-trivially affected" << endl;
                        // TODO
                    }
                }
            }
        }
        assert(affected_abstractions.size() == 1);
        atomic_symmetries_by_affected_abs[*affected_abstractions.begin()].push_back(gen_index);
    }
    return true;
}

bool Symmetries::find_to_be_merged_abstractions(const vector<Abstraction *>& abstractions,
                                                set<int> &abs_to_merge) {
    assert(abs_to_merge.empty());
    vector<set<int> > non_trivially_affected_abstractions;
    vector<bool> atomic_symmetries;
    vector<vector<int> > atomic_symmetries_by_affected_abs;
    bool found_symmetry = find_symmetries(abstractions, non_trivially_affected_abstractions, atomic_symmetries, atomic_symmetries_by_affected_abs, false);
    if (found_symmetry) {
        int smallest_symmetry_index = -1;
        int smallest_symmetry_size = numeric_limits<int>::max();
        for (size_t i = 0; i < non_trivially_affected_abstractions.size(); ++i) {
            // find the smallest symmetry (with the least number of affected abstractions)
            int number_affected_abs = non_trivially_affected_abstractions[i].size();
            assert(number_affected_abs > 0);
            if (number_affected_abs < smallest_symmetry_size) {
                smallest_symmetry_size = number_affected_abs;
                smallest_symmetry_index = i;
            }
        }
        assert(smallest_symmetry_index != -1);
        abs_to_merge.insert(non_trivially_affected_abstractions[smallest_symmetry_index].begin(),
                            non_trivially_affected_abstractions[smallest_symmetry_index].end());
    }
    return found_symmetry;
}

bool Symmetries::find_symmetries(const vector<Abstraction *>& abstractions,
                                 vector<set<int> > &non_trivially_affected_abstractions,
                                 vector<bool> &atomic_symmetries,
                                 vector<vector<int> > &atomic_symmetries_by_affected_abs,
                                 bool find_atomic_symmetry) {
    /**
     * Find symmetries for abstractions. If atomic is true, find abstraction stabilized
     * symmetries and check wether they are atomic or not. If atomic is false, find
     * general symmetries. In both cases, non_trivially_affected_abstractions contain
     * the abstraction indices that are affected by the symmetries. Returns true if any
     * symmetry is found at all or false.
     */
    assert(non_trivially_affected_abstractions.empty());
    cout << "Computing generators for " << (find_atomic_symmetry? "" : "non ") << "stabilized symmetries" << endl;
    gc.compute_generators(abstractions, find_atomic_symmetry);

    if (get_num_generators() == 0) {
        cout << "No generators found! Done searching for symmetries. [t=" << g_timer << "]" << endl;
        return false;
    }

    int num_abstractions = get_permutations_wrapper().num_abstractions;
    unsigned int num_generators = get_num_generators();

    // Compute all non-trivially affected abstractions:
    non_trivially_affected_abstractions.resize(num_generators, set<int>());
    atomic_symmetries.resize(num_generators, false);
    atomic_symmetries_by_affected_abs.resize(abstractions.size(), vector<int>());
    bool found_atomic_symmetry = false;
    for (int gen_index = 0; gen_index < num_generators; ++gen_index) {
        cout << "generator " << gen_index << endl;
        // Find all abstractions not mapped to themselves
        set<int> &affected_abs = non_trivially_affected_abstractions[gen_index];
        bool atomic_symmetry = true;
        for (unsigned int index = 0; index < num_abstractions; ++index) {
            if (abstractions[index] == 0)
                continue;
            unsigned int to_index = get_generator(gen_index)->get_value(index);
            if (index == to_index)
                continue;
            if (atomic_symmetry) {
                // index != to_index, thus the symmetry is not atomic
                atomic_symmetry = false;
                if (find_atomic_symmetry) {
                    // if only interested in atomic symmetries, stop here
                    break;
                }
            }
            affected_abs.insert(index);
            cout << "abstraction " << abstractions[index]->description() << " mapped to " << abstractions[to_index]->description() << endl;
        }
        if (find_atomic_symmetry && !atomic_symmetry) {
            // if only interested in atomic symmetries and symmetry is not atomic, skip symmetry
            cout << "symmetry " << gen_index << "not atomic, skipping" << endl;
            // atomic_symmetries contains a default false already
            continue;
        }

        // Find all abstractions whose states are not all mapped to states from the same abstraction
        for (unsigned int index = num_abstractions; index < get_permutations_wrapper().num_abs_and_states; ++index) {
            int abs_index = get_permutations_wrapper().get_var_by_index(index);
            if (abstractions[abs_index] == 0)
                continue;
            unsigned int to_index = get_generator(gen_index)->get_value(index);
            if (index == to_index)
                continue;
            if (affected_abs.insert(abs_index).second) {
                cout << "abstraction " << abstractions[abs_index]->description() << " non-trivially affected" << endl;
                if (atomic_symmetry && affected_abs.size() > 1) {
                    // more than one abstraction is affected
                    atomic_symmetry = false;
                    if (find_atomic_symmetry) {
                        // if only interested in atomic symmetries, stop here
                        break;
                    }
                }
            }
        }
        cout << "It is " << (atomic_symmetry ? "" : "not ") << "an atomic symmetry" << endl;
        if (find_atomic_symmetry && atomic_symmetry)
            found_atomic_symmetry = true;
        if (atomic_symmetry) {
            assert(affected_abs.size() == 1);
            atomic_symmetries_by_affected_abs[*affected_abs.begin()].push_back(gen_index);
//            for (unsigned int index = 0; index < get_permutations_wrapper().num_abs_and_states; ++index) {
//                int abs_index = get_permutations_wrapper().get_var_by_index(index);
//                if (abstractions[abs_index] == 0) {
//                    cout << "abstraction not in use" << endl;
//                    cout << "but would do: " << endl;
//                    //continue;
//                }
//                unsigned int to_index = get_generator(gen_index)->get_value(index);
//                cout << index << " mapped to " << to_index << endl;
//            }
        }
        atomic_symmetries[gen_index] = atomic_symmetry;
    }
    if (find_atomic_symmetry)
        return found_atomic_symmetry;
    // found symmetries
    return true;
}

void Symmetries::apply_symmetry(const vector<Abstraction *> &abstractions, int generator_index) const {
    if (get_num_generators() == 0) {
        cerr << "You first have to find symmetries before you can apply one of them!" << endl;
        exit_with(EXIT_CRITICAL_ERROR);
    }
    assert(0 <= generator_index && generator_index < get_num_generators());
    cout << "Creating equivalence relations from symmetries. [t=" << g_timer << "]" << endl;

    // Abstracting by the generators as follows:
    // Creating a graph with nodes being the abstract states and the edges represent the connections as given by the generators.
    // It seems like this process should better be done in all abstractions in parallel,
    // since we could exploit all the compositions of these abstractions this way as well.
    // Later we can reduce the first part - the abstraction vars and save some space/time
    unsigned int num_states = get_permutations_wrapper().num_abs_and_states;
    unsigned int num_abstractions = get_permutations_wrapper().num_abstractions;
    // The graph is represented by vector of to_nodes for each node. (Change to sets?)
    vector<vector<unsigned int> > graph(num_states, vector<unsigned int>());
    for (unsigned int index = num_abstractions; index < get_permutations_wrapper().num_abs_and_states; ++index) {
        int abs_index = get_permutations_wrapper().get_var_by_index(index);
        if (abstractions[abs_index] == 0)
            continue;
        unsigned int to_index = get_generator(generator_index)->get_value(index);
        if (index == to_index)
            continue;
        graph[index].push_back(to_index);
    }
    SCC scc(graph);
    const vector<vector<unsigned int> >& result = scc.get_result();

    // TODO: check if we really need this complex (at least it seems to be complex)
    // method of generating equivalence classes in the case of only applying one generator
    // Generate final result. Going over the result, putting the nodes to their respective places.
    vector<EquivRel> equivalence_relations(abstractions.size(), EquivRel());
    for (size_t i = 0; i < abstractions.size(); ++i) {
        if (abstractions[i] == 0)  //In case the abstraction is empty
            continue;

        equivalence_relations[i].resize(result.size());
    }
    for (unsigned int eqiv=0; eqiv < result.size(); eqiv++) {
        for (unsigned int i=0; i < result[eqiv].size(); i++) {
            unsigned int idx = result[eqiv][i];
            if (idx < get_permutations_wrapper().num_abstractions)
                continue;
            pair<int, AbstractStateRef> vals = get_permutations_wrapper().get_var_val_by_index(idx);
            equivalence_relations[vals.first][eqiv].push_front(vals.second);
        }
    }
    // Then, going over the outcome, removing empty equivalence classes.
    for (unsigned int eqiv=result.size(); eqiv > 0; --eqiv) {
        for (size_t i = 0; i < abstractions.size(); ++i) {
            if (abstractions[i] == 0)  //In case the abstraction is empty
                continue;

            if (equivalence_relations[i][eqiv - 1].size() == 0)
                equivalence_relations[i].erase(equivalence_relations[i].begin() + eqiv - 1);
        }
    }

    cout << "==========================================================================================" << endl;
    cout << "Abstracting everything by the equivalence relations. [t=" << g_timer << "]" << endl;

    for (size_t i = 0; i < abstractions.size(); ++i) {
        if (abstractions[i] == 0)  //In case the abstraction is empty
            continue;

        // Abstracting by the computed equivalence relation
        if (equivalence_relations[i].size() > abstractions[i]->size()) {
            cerr << "Something is seriously wrong here!!" << endl;
            exit_with(EXIT_CRITICAL_ERROR);
        }
        if (equivalence_relations[i].size() == abstractions[i]->size()) {
            cout << abstractions[i]->tag() << " not abstracted due to symmetries." << endl;
            continue;
        }
//      cout << "Abstracting from " << abstractions[i]->size() << " to " << equivalence_relations[i].size() << " states!" << endl;
        abstractions[i]->apply_abstraction(equivalence_relations[i]);
    }
    cout << "Done abstracting. [t=" << g_timer << "]" << endl;
    cout << "==========================================================================================" << endl;
}

void Symmetries::apply_symmetries(const vector<Abstraction *> &abstractions, const vector<int> &indices) const {
    if (get_num_generators() == 0) {
        cerr << "You first have to find symmetries before you can apply one of them!" << endl;
        exit_with(EXIT_CRITICAL_ERROR);
    }
    cout << "Creating equivalence relations from symmetries. [t=" << g_timer << "]" << endl;

    // Abstracting by the generators as follows:
    // Creating a graph with nodes being the abstract states and the edges represent the connections as given by the generators.
    // It seems like this process should better be done in all abstractions in parallel,
    // since we could exploit all the compositions of these abstractions this way as well.
    // Later we can reduce the first part - the abstraction vars and save some space/time
    unsigned int num_states = get_permutations_wrapper().num_abs_and_states;
    unsigned int num_abstractions = get_permutations_wrapper().num_abstractions;
    // The graph is represented by vector of to_nodes for each node. (Change to sets?)
    vector<vector<unsigned int> > graph(num_states, vector<unsigned int>());
    for (unsigned int index = num_abstractions; index < get_permutations_wrapper().num_abs_and_states; ++index) {
        int abs_index = get_permutations_wrapper().get_var_by_index(index);
        if (abstractions[abs_index] == 0)
            continue;
        for (int i = 0; i < indices.size(); ++i) {
            // Going over the generators, for each just add the edges.
            unsigned int to_index = get_generator(indices[i])->get_value(index);
            if (index == to_index)
                continue;
            graph[index].push_back(to_index);
        }
    }
    SCC scc(graph);
    const vector<vector<unsigned int> >& result = scc.get_result();

    // TODO: check if we really need this complex (at least it seems to be complex)
    // method of generating equivalence classes in the case of only applying one generator
    // Generate final result. Going over the result, putting the nodes to their respective places.
    vector<EquivRel> equivalence_relations(abstractions.size(), EquivRel());
    for (size_t i = 0; i < abstractions.size(); ++i) {
        if (abstractions[i] == 0)  //In case the abstraction is empty
            continue;

        equivalence_relations[i].resize(result.size());
    }
    for (unsigned int eqiv=0; eqiv < result.size(); eqiv++) {
        for (unsigned int i=0; i < result[eqiv].size(); i++) {
            unsigned int idx = result[eqiv][i];
            if (idx < get_permutations_wrapper().num_abstractions)
                continue;
            pair<int, AbstractStateRef> vals = get_permutations_wrapper().get_var_val_by_index(idx);
            equivalence_relations[vals.first][eqiv].push_front(vals.second);
        }
    }
    // Then, going over the outcome, removing empty equivalence classes.
    for (unsigned int eqiv=result.size(); eqiv > 0; --eqiv) {
        for (size_t i = 0; i < abstractions.size(); ++i) {
            if (abstractions[i] == 0)  //In case the abstraction is empty
                continue;

            if (equivalence_relations[i][eqiv - 1].size() == 0)
                equivalence_relations[i].erase(equivalence_relations[i].begin() + eqiv - 1);
        }
    }

    cout << "==========================================================================================" << endl;
    cout << "Abstracting everything by the equivalence relations. [t=" << g_timer << "]" << endl;

    for (size_t i = 0; i < abstractions.size(); ++i) {
        if (abstractions[i] == 0)  //In case the abstraction is empty
            continue;

        // Abstracting by the computed equivalence relation
        if (equivalence_relations[i].size() > abstractions[i]->size()) {
            cerr << "Something is seriously wrong here!!" << endl;
            exit_with(EXIT_CRITICAL_ERROR);
        }
        if (equivalence_relations[i].size() == abstractions[i]->size()) {
            cout << abstractions[i]->tag() << " not abstracted due to symmetries." << endl;
            continue;
        }
//      cout << "Abstracting from " << abstractions[i]->size() << " to " << equivalence_relations[i].size() << " states!" << endl;
        abstractions[i]->apply_abstraction(equivalence_relations[i]);
    }
    cout << "Done abstracting. [t=" << g_timer << "]" << endl;
    cout << "==========================================================================================" << endl;
}

void Symmetries::print_generators_stat() const {
    cout << "Found generators:" << endl;
    int cycles_1 = 0, cycles_2 = 0, cycles_3 = 0, cycles_more = 0;
    int non_linear_merges_1 = 0, non_linear_merges_2 = 0, non_linear_merges_3 = 0, non_linear_merges_more = 0;
    int linear_merges_none = 0, linear_merges_2 = 0, linear_merges_3 = 0, linear_merges_more = 0;
    for (size_t i = 0; i < get_num_generators(); ++i) {
        const Permutation* perm = get_generator(i);
        perm->print_cycle_notation();
        perm->print_variables_by_cycles();
        cout << "Generator " << i << " has maximal variable cycle of size " << perm->get_maximal_variables_cycle_size() << endl;
        if (perm->get_maximal_variables_cycle_size() < 2) {
            cycles_1++;
        } else if (perm->get_maximal_variables_cycle_size() == 2) {
            cycles_2++;
            // Getting merge statistics
            int lin_var_merge = perm->calculate_number_variables_to_merge(true);
            cout << "Generator " << i << " for linear merge demands " << lin_var_merge << " merges" << endl;
            int non_lin_var_merge = perm->calculate_number_variables_to_merge(false);
            cout << "Generator " << i << " for non linear merge demands " << non_lin_var_merge << " merges" << endl;

            if (lin_var_merge < 2) {
                linear_merges_none++;
            } else if (lin_var_merge == 2) {
                linear_merges_2++;
            } else if (lin_var_merge == 3) {
                linear_merges_3++;
            } else if (lin_var_merge > 3) {
                linear_merges_more++;
            }

            if (non_lin_var_merge == 1) {
                non_linear_merges_1++;
            } else if (non_lin_var_merge == 2) {
                non_linear_merges_2++;
            } else if (non_lin_var_merge == 3) {
                non_linear_merges_3++;
            } else if (non_lin_var_merge > 3) {
                non_linear_merges_more++;
            }

        } else if (perm->get_maximal_variables_cycle_size() == 3) {
            cycles_3++;
        } else if (perm->get_maximal_variables_cycle_size() > 3) {
            cycles_more++;
        }
    }
    cout << "Number of generators with no cycles is " << cycles_1 << endl;
    cout << "Number of generators with maximal cycle size 2 is " << cycles_2 << endl;
    cout << "Number of generators with maximal cycle size 3 is " << cycles_3 << endl;
    cout << "Number of generators with maximal cycle size greater than 3 is " << cycles_more << endl;

    cout << "Number of generators that do not demand any merge is " << linear_merges_none << endl;
    cout << "Number of generators that demand linear merge of size 2 is " << linear_merges_2 << endl;
    cout << "Number of generators that demand linear merge of size 3 is " << linear_merges_3 << endl;
    cout << "Number of generators that demand linear merge of size greater than 3 is " << linear_merges_more << endl;

    cout << "Number of generators that demand non linear merge of size 1 is " << non_linear_merges_1 << endl;
    cout << "Number of generators that demand non linear merge of size 2 is " << non_linear_merges_2 << endl;
    cout << "Number of generators that demand non linear merge of size 3 is " << non_linear_merges_3 << endl;
    cout << "Number of generators that demand non linear merge of size greater than 3 is " << non_linear_merges_more << endl;

}

const Permutation* Symmetries::get_generator(int ind) const {
    assert(ind >= 0 && ind < get_num_generators());
    return gc.get_generators()[ind];
}
