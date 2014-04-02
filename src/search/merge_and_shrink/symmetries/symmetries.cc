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

bool Symmetries::is_atomar_generator(const vector<Abstraction *> abstractions, int gen_index) const {
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

bool Symmetries::find_and_apply_atomar_symmetries(const vector<Abstraction *> &abstractions) {
    cout << "Applying all atomar symmetries:" << endl;
    bool applied_symmetry = false;
    while (true) {
        vector<set<int> > non_trivially_affected_abstractions;
        vector<bool> atomar_symmetries;
        vector<vector<int> > atomar_symmetries_by_affected_abs;
        bool found_symmetry = find_symmetries(abstractions, non_trivially_affected_abstractions, atomar_symmetries, atomar_symmetries_by_affected_abs, true);
        if (!found_symmetry) {
            if (applied_symmetry)
                cout << "Applied at least one atomar symmetry." << endl;
            else
                cout << "Did not find any applicable atomar symmetries." << endl;
            return applied_symmetry;
        }
        if (!applied_symmetry) {
            applied_symmetry = true;
        }
        int num_atomar_sym = 0;
        for (size_t i = 0; i < atomar_symmetries.size(); ++i) {
            if (atomar_symmetries[i])
                ++num_atomar_sym;
        }
        if (version == 1) {
            int most_affected_abstraction_index = -1;
            int most_affected_abstraction_size = 0;
            for (size_t abs_index = 0; abs_index < abstractions.size(); ++abs_index) {
                if (abstractions[abs_index] == 0) {
                    continue;
                }
                int no_of_atomar_symmetries_affected = atomar_symmetries_by_affected_abs[abs_index].size();
                if (no_of_atomar_symmetries_affected > most_affected_abstraction_size) {
                    most_affected_abstraction_size = no_of_atomar_symmetries_affected;
                    most_affected_abstraction_index = abs_index;
                }
            }
            assert(most_affected_abstraction_index != -1);
            assert(most_affected_abstraction_size > 0);
            // TODO: apply all atomic abstractions combined rather than just the "biggest one"?
            // See conjecture below.
            apply_symmetries(abstractions, atomar_symmetries_by_affected_abs[most_affected_abstraction_index]);
        } else {
            for (size_t gen_index = 0; gen_index < atomar_symmetries.size(); ++gen_index) {
                if (!atomar_symmetries[gen_index]) {
                    continue;
                }
                // TODO: is it correct that several atomar symmetries can always be combined?
                // Conjecture: the resulting atomar symmetry would not be atomar in our
                // definition's sense, but the problematic case why we need atomar
                // symmetries and cannot use the more general "local" symmetries cannot
                // arise.
                bool atomar_generator = is_atomar_generator(abstractions, gen_index);
                assert(atomar_generator);
                apply_symmetry(abstractions, gen_index);
//                cout << "SAFETY CHECK..." << endl;
//                cout << "original number of generators: " << num_atomar_sym << endl;
//                cout << "just applied " << gen_index << endl;
//                Symmetries symmetries(*this);
//                vector<set<int> > non_trivially_affected_abstractions_other;
//                vector<bool> atomar_symmetries_other;
//                symmetries.find_symmetries(abstractions, non_trivially_affected_abstractions_other, atomar_symmetries_other, true);
//                int num_atomar_sym_other = 0;
//                for (size_t i = 0; i < atomar_symmetries_other.size(); ++i) {
//                    if (atomar_symmetries_other[i])
//                        ++num_atomar_sym_other;
//                }
//                int num_remaining_sym = num_atomar_sym - (gen_index + 1);
//                cout << "found atomar symmetries: " << num_atomar_sym_other << endl;
//                cout << "remaining symmetries: " << num_remaining_sym << endl;
//                //assert(num_atomar_sym_other == num_remaining_sym);
//                cout << "now we have " << num_atomar_sym_other << " remaining symmetries" << endl;
//                for (size_t i = 0; i < non_trivially_affected_abstractions_other.size(); ++i) {
//                    for (set<int>::iterator it = non_trivially_affected_abstractions_other[i].begin();
//                         i != non_trivially_affected_abstractions_other[i].end(); ++it) {
//                    }
//                }
//                cout << "DONE" << endl;
            }
        }
    }
}

bool Symmetries::find_to_be_merged_abstractions(const vector<Abstraction *>& abstractions,
                                                set<int> &abs_to_merge) {
    assert(abs_to_merge.empty());
    vector<set<int> > non_trivially_affected_abstractions;
    vector<bool> atomar_symmetries;
    vector<vector<int> > atomar_symmetries_by_affected_abs;
    bool found_symmetry = find_symmetries(abstractions, non_trivially_affected_abstractions, atomar_symmetries, atomar_symmetries_by_affected_abs, false);
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
                                 vector<bool> &atomar_symmetries,
                                 vector<vector<int> > &atomar_symmetries_by_affected_abs,
                                 bool find_atomar_symmetry) {
    /**
     * Find symmetries for abstractions. If atomar is true, find abstraction stabilized
     * symmetries and check wether they are atomar or not. If atomar is false, find
     * general symmetries. In both cases, non_trivially_affected_abstractions contain
     * the abstraction indices that are affected by the symmetries. Returns true if any
     * symmetry is found at all or false.
     */
    assert(non_trivially_affected_abstractions.empty());
    cout << "Computing generators for " << (find_atomar_symmetry? "" : "non ") << "stabilized symmetries" << endl;
    gc.compute_generators(abstractions, find_atomar_symmetry);

    if (get_num_generators() == 0) {
        cout << "No generators found! Done searching for symmetries. [t=" << g_timer << "]" << endl;
        return false;
    }

    unsigned int num_abstractions = get_permutations_wrapper().num_abstractions;
    unsigned int num_generators = get_num_generators();

    // Compute all non-trivially affected abstractions:
    non_trivially_affected_abstractions.resize(num_generators, set<int>());
    atomar_symmetries.resize(num_generators, false);
    atomar_symmetries_by_affected_abs.resize(abstractions.size(), vector<int>());
    bool found_atomar_symmetry = false;
    for (int gen_index = 0; gen_index < num_generators; ++gen_index) {
        cout << "generator " << gen_index << endl;
        // Find all abstractions not mapped to themselves
        vector<vector<unsigned int> > graph(num_abstractions, vector<unsigned int>());
        set<int> &affected_abs = non_trivially_affected_abstractions[gen_index];
        bool atomar_symmetry = true;
        for (unsigned int index = 0; index < num_abstractions; ++index) {
            if (abstractions[index] == 0)
                continue;
            unsigned int to_index = get_generator(gen_index)->get_value(index);
            if (index == to_index)
                continue;
            if (atomar_symmetry) {
                // index != to_index, thus the symmetry is not atomar
                atomar_symmetry = false;
                if (find_atomar_symmetry) {
                    // if only interested in atomar symmetries, stop here
                    break;
                }
            }
            affected_abs.insert(index);
            cout << "abstraction " << abstractions[index]->description() << " mapped to " << abstractions[to_index]->description() << endl;
            graph[index].push_back(to_index);
        }
        if (find_atomar_symmetry && !atomar_symmetry) {
            // if only interested in atomar symmetries and symmetry is not atomar, skip symmetry
            cout << "symmetry " << gen_index << "not atomar, skipping" << endl;
            // atomar_symmetries contains a default false already
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
                if (atomar_symmetry && affected_abs.size() > 1) {
                    // more than one abstraction is affected
                    atomar_symmetry = false;
                    if (find_atomar_symmetry) {
                        // if only interested in atomar symmetries, stop here
                        break;
                    }
                }
            }
        }
        cout << "It is " << (atomar_symmetry ? "" : "not ") << "an atomar symmetry" << endl;
        if (find_atomar_symmetry && atomar_symmetry)
            found_atomar_symmetry = true;
        if (atomar_symmetry) {
            assert(affected_abs.size() == 1);
            atomar_symmetries_by_affected_abs[*affected_abs.begin()].push_back(gen_index);
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
        atomar_symmetries[gen_index] = atomar_symmetry;
    }
    if (find_atomar_symmetry)
        return found_atomar_symmetry;
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
        if (!abstractions[i]->are_distances_computed())
            abstractions[i]->compute_distances();
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
        if (!abstractions[i]->are_distances_computed())
            abstractions[i]->compute_distances();
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
