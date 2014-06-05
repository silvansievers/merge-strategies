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
    : gc(options),
      type_of_symmetries(TypeOfSymmetries(options.get_enum("type_of_symmetries"))),
      build_stabilized_pdg(options.get<bool>("build_stabilized_pdg")),
      atomic_symmetries(0),
      binary_symmetries(0),
      other_symmetries(0) {
}

int Symmetries::find_and_apply_symmetries(const vector<Abstraction *>& abstractions,
                                          set<int> &abs_to_merge) {
    assert(abs_to_merge.empty());
    int number_of_applied_symmetries = 0;
    while (true) {
        vector<set<int> > affected_abstractions_by_generator;
        vector<set<int> > mapped_abstractions_by_generator;
        vector<int> atomic_generators;
        vector<int> local_generators;
        bool found_symmetry = find_symmetries(abstractions,
                                              affected_abstractions_by_generator,
                                              mapped_abstractions_by_generator,
                                              atomic_generators,
                                              local_generators);
        if (found_symmetry) {
            for (size_t generator = 0; generator < affected_abstractions_by_generator.size(); ++generator) {
                switch (affected_abstractions_by_generator[generator].size()) {
                case 0:
                    cerr << "Found an identity generator!" << endl;
                    exit_with(EXIT_CRITICAL_ERROR);
                    break;
                case 1:
                    ++atomic_symmetries;
                    break;
                case 2:
                    ++binary_symmetries;
                    break;
                default:
                    ++other_symmetries;
                    break;
                }
            }

            switch (type_of_symmetries) {
            case ATOMIC:
                if (atomic_generators.empty()) {
                    // find the smallest symmetry (with the least number of affected abstractions)
                    find_smallest_generator(affected_abstractions_by_generator, abs_to_merge);
                    return number_of_applied_symmetries;
                } else {
                    ++number_of_applied_symmetries;
                    apply_symmetries(abstractions, atomic_generators);
                }
                break;
            case LOCAL:
                if (local_generators.empty()) {
                    if (!mapped_abstractions_by_generator.empty()) {
                        // find the smallest symmetry (with the least number of mapped abstractions)
                        find_smallest_generator(mapped_abstractions_by_generator, abs_to_merge);
                    } else {
                        // find the smallest symmetry (with the least number of affected abstractions)
                        find_smallest_generator(affected_abstractions_by_generator, abs_to_merge);
                    }
                    return number_of_applied_symmetries;
                } else {
                    ++number_of_applied_symmetries;
                    apply_symmetries(abstractions, local_generators);
                }
                break;
            case ANY:
                ++number_of_applied_symmetries;
                if (get_num_generators() == 0) {
                    return number_of_applied_symmetries;
                } else {
                    ++number_of_applied_symmetries;
                    vector<int> all_generators;
                    all_generators.reserve(get_num_generators());
                    for (size_t i = 0; i < get_num_generators(); ++i)
                        all_generators.push_back(i);
                    apply_symmetries(abstractions, all_generators);
                }
                break;
            }
        } else {
            break;
        }
    }
    return number_of_applied_symmetries;
}

void Symmetries::find_smallest_generator(const vector<set<int> > &abstractions_by_generator,
                                         set<int> &abs_to_merge) const {
    int smallest_generator_index = -1;
    int smallest_generator_size = numeric_limits<int>::max();
    for (size_t generator = 0; generator < abstractions_by_generator.size(); ++generator) {
        int number_affected_abs = abstractions_by_generator[generator].size();
        assert(number_affected_abs > 0);
        if (number_affected_abs < smallest_generator_size) {
            smallest_generator_size = number_affected_abs;
            smallest_generator_index = generator;
        }
    }
    assert(smallest_generator_index != -1);
    abs_to_merge.insert(abstractions_by_generator[smallest_generator_index].begin(),
                        abstractions_by_generator[smallest_generator_index].end());
}

bool Symmetries::find_symmetries(const vector<Abstraction *>& abstractions,
                                 vector<set<int> > &affected_abstractions_by_generator,
                                 vector<set<int> > &mapped_abstractions_by_generator,
                                 vector<int> &atomic_generators,
                                 vector<int> &local_generators) {
    // Find (non) abstraction stabilized symmetries for abstractions depending
    // on the chosen option.
    // When returning, affected_abstractions_by_generator contains the
    // abstraction indices that are affected by each generator,
    // atomic_generators contains the indices of atomic generators, and
    // local_generators contains the indices of local generators.
    // Returns true if any symmetry is found at all and false otherwise.

    assert(affected_abstractions_by_generator.empty());
    assert(mapped_abstractions_by_generator.empty());

    // We must make sure that all abstractions distances have been computed
    // because of the nasty possible side effect of pruning irrelevant
    // states and because the application of an equivalence realtion to an
    // abstraction requires distances to be computed.
    for (size_t i = 0; i < abstractions.size(); ++i) {
        if (abstractions[i])
            abstractions[i]->compute_distances();
    }

    gc.compute_generators(abstractions);

    unsigned int num_generators = get_num_generators();
    if (num_generators == 0) {
        cout << "No generators found! Done searching for symmetries. [t=" << g_timer << "]" << endl;
        return false;
    }

    int num_abstractions = get_pw().num_abstractions;

    // Find all affected abstractions
    affected_abstractions_by_generator.resize(num_generators, set<int>());
    mapped_abstractions_by_generator.resize(num_generators, set<int>());
    for (int gen_index = 0; gen_index < num_generators; ++gen_index) {
        cout << "generator " << gen_index << endl;
        set<int> &affected_abs = affected_abstractions_by_generator[gen_index];
        set<int> &mapped_abs = mapped_abstractions_by_generator[gen_index];
        bool local_generator = true;

        if (!build_stabilized_pdg) {
            // Find all abstractions not mapped to themselves
            for (unsigned int index = 0; index < num_abstractions; ++index) {
                if (abstractions[index]) {
                    unsigned int to_index = get_generator(gen_index)->get_value(index);
                    if (index != to_index) {
                        local_generator = false;
                        affected_abs.insert(index);
                        // we could potentially also add to_index here, but
                        // it will must be added when considering to_index
                        // anyway.
                        mapped_abs.insert(index);
                        cout << "abstraction " << abstractions[index]->description()
                             << " mapped to " << abstractions[to_index]->description() << endl;
                    }
                }
            }
        }

        // Find all abstractions whose states are not mapped to themselves.
        for (unsigned int index = num_abstractions; index < get_pw().num_abs_and_states; ++index) {
            int abs_index = get_pw().get_var_by_index(index);
            if (!abstractions[abs_index]) {
                cerr << "found an abstract state belonging to an invalid abstraction" << endl;
                exit_with(EXIT_CRITICAL_ERROR);
            }
            unsigned int to_index = get_generator(gen_index)->get_value(index);
            if (index != to_index) {
                if (affected_abs.insert(abs_index).second) {
                    cout << "abstraction " << abstractions[abs_index]->description() << " is affected" << endl;
                }
            }
        }

        assert(mapped_abs.size() <= affected_abs.size());
        if (affected_abs.size() == 1) {
            assert(mapped_abs.empty());
            atomic_generators.push_back(gen_index);
        } else if (affected_abs.size() < 1) {
            cerr << "Something is wrong! The generator is the identity generator." << endl;
            exit_with(EXIT_CRITICAL_ERROR);
        }
        if (local_generator) {
            local_generators.push_back(gen_index);
        }
    }
    return true;
}

void Symmetries::apply_symmetries(const vector<Abstraction *> &abstractions,
                                  const vector<int> &generator_indices) const {
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
    unsigned int num_states = get_pw().num_abs_and_states;
    unsigned int num_abstractions = get_pw().num_abstractions;
    // The graph is represented by vector of to_nodes for each node. (Change to sets?)
    vector<vector<unsigned int> > graph(num_states, vector<unsigned int>());
    for (unsigned int index = num_abstractions; index < get_pw().num_abs_and_states; ++index) {
        int abs_index = get_pw().get_var_by_index(index);
        if (abstractions[abs_index]) {
            for (size_t i = 0; i < generator_indices.size(); ++i) {
                // Going over the generators, for each just add the edges.
                unsigned int to_index = get_generator(generator_indices[i])->get_value(index);
                if (index != to_index)
                    graph[index].push_back(to_index);
            }
        }
    }
    SCC scc(graph);
    const vector<vector<unsigned int> >& result = scc.get_result();

    // TODO: check if we really need this complex (at least it seems to be complex)
    // method of generating equivalence classes in the case of only applying one generator
    // Generate final result. Going over the result, putting the nodes to their respective places.
    vector<EquivRel> equivalence_relations(abstractions.size(), EquivRel());
    for (size_t abs_index = 0; abs_index < abstractions.size(); ++abs_index) {
        if (abstractions[abs_index])
            equivalence_relations[abs_index].resize(result.size());
    }
    for (unsigned int eqiv=0; eqiv < result.size(); eqiv++) {
        for (unsigned int i=0; i < result[eqiv].size(); i++) {
            unsigned int idx = result[eqiv][i];
            if (idx < get_pw().num_abstractions)
                continue;
            pair<int, AbstractStateRef> vals = get_pw().get_var_val_by_index(idx);
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

const Permutation* Symmetries::get_generator(int ind) const {
    assert(ind >= 0 && ind < get_num_generators());
    return gc.get_generators()[ind];
}
