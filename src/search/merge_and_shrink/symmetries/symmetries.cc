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

pair<int, int> Symmetries::find_and_apply_symmetries(vector<Abstraction *> &abstractions,
                                                     set<int> &abs_to_merge) {
    assert(abs_to_merge.empty());
    int number_of_applied_symmetries = 0;
    int number_of_collapsed_abstractions = 0;
    while (true) {
        //for (size_t i = 0; i < abstractions.size(); ++i) {
        //    cout << i << " " << abstractions[i] << endl;
        //}

        vector<set<int> > affected_abstractions_by_generator;
        vector<set<int> > mapped_abstractions_by_generator;
        vector<set<int> > affected_not_mapped_abs_by_gen;
        vector<vector<vector<int> > > cycles_by_generators;
        vector<int> atomic_generators;
        vector<int> local_generators;
        bool found_symmetry = find_symmetries(abstractions,
                                              affected_abstractions_by_generator,
                                              mapped_abstractions_by_generator,
                                              affected_not_mapped_abs_by_gen,
                                              cycles_by_generators,
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
                    return make_pair(number_of_applied_symmetries, 0);
                } else {
                    ++number_of_applied_symmetries;
                    apply_symmetries(abstractions, atomic_generators);
                }
                break;
            case LOCAL:
                if (local_generators.empty()) {
                    bool mapped_abstractions = false;
                    for (size_t i = 0; i < mapped_abstractions_by_generator.size(); ++i) {
                        if (!mapped_abstractions_by_generator[i].empty()) {
                            mapped_abstractions = true;
                            break;
                        }
                    }
                    if (mapped_abstractions) {
                        // find the smallest symmetry (with the least number of mapped abstractions)
                        find_smallest_generator(mapped_abstractions_by_generator, abs_to_merge);
                    } else {
                        // find the smallest symmetry (with the least number of affected abstractions)
                        find_smallest_generator(affected_abstractions_by_generator, abs_to_merge);
                    }
                    return make_pair(number_of_applied_symmetries, 0);
                } else {
                    ++number_of_applied_symmetries;
                    apply_symmetries(abstractions, local_generators);
                }
                break;
            case AT_MOST_ONE_CYCLE:
                if (get_num_generators() == 0) {
                    return make_pair(number_of_applied_symmetries, number_of_collapsed_abstractions);
                } else {
                    int largest_cycle_generator_index = -1;
                    int largest_cycle_index = -1;
                    int largest_cycle_size = 0;
                    for (size_t generator = 0; generator < cycles_by_generators.size(); ++generator) {
                        const vector<vector<int> > &cycles = cycles_by_generators[generator];
                        if (cycles.size() == 1 &&
                                affected_not_mapped_abs_by_gen[generator].empty()) {
                            // only consider generators that have only one cycle
                            // of mapped abstractions and do not affect others.
                            //if (!cycles.empty())
                            //    cout << "generator " << generator << endl;
                            for (size_t i = 0; i < cycles.size(); ++i) {
                                const vector<int> &cycle = cycles[i];
                                /*cout << "cycle " << i << endl;
                                for (size_t j = 0; j < cycle.size(); ++j) {
                                    cout << cycle[j] << " ";
                                }
                                cout << endl;*/

                                if (cycle.size() > largest_cycle_size) {
                                    largest_cycle_size = cycle.size();
                                    largest_cycle_index = i;
                                    largest_cycle_generator_index = generator;
                                }
                            }
                        }
                    }
                    if (largest_cycle_generator_index != -1) {
                        ++number_of_applied_symmetries;
                        assert(largest_cycle_index == 0);
                        const vector<int> &collapsed_abs =
                                cycles_by_generators[largest_cycle_generator_index][largest_cycle_index];
                        Abstraction *chosen_representative = abstractions[collapsed_abs[0]];
                        for (size_t abs_index = 1; abs_index < collapsed_abs.size(); ++abs_index) {
                            Abstraction *abs = abstractions[abs_index];
                            chosen_representative->merge_abstraction_into(abs);
                            chosen_representative->normalize();
                            abs->release_memory();
                            abstractions[abs_index] = 0;
                        }
                        number_of_collapsed_abstractions += collapsed_abs.size() - 1;
                    } else {
                        bool mapped_abstractions = false;
                        for (size_t i = 0; i < mapped_abstractions_by_generator.size(); ++i) {
                            if (!mapped_abstractions_by_generator[i].empty()) {
                                mapped_abstractions = true;
                                break;
                            }
                        }
                        if (mapped_abstractions) {
                            // find the smallest symmetry (with the least number of mapped abstractions)
                            find_smallest_generator(mapped_abstractions_by_generator, abs_to_merge);
                        } else {
                            // find the smallest symmetry (with the least number of affected abstractions)
                            find_smallest_generator(affected_abstractions_by_generator, abs_to_merge);
                        }
                        return make_pair(number_of_applied_symmetries, number_of_collapsed_abstractions);
                    }
                }
                break;
            case ANY:
                cerr << "not implemented" << endl;
                exit_with(EXIT_CRITICAL_ERROR);
                return make_pair(number_of_applied_symmetries, 0);
                break;
            }
        } else {
            break;
        }
    }
    return make_pair(number_of_applied_symmetries, number_of_collapsed_abstractions);
}

void Symmetries::find_smallest_generator(const vector<set<int> > &abstractions_by_generator,
                                         set<int> &abs_to_merge) const {
    int smallest_generator_index = -1;
    int smallest_generator_size = numeric_limits<int>::max();
    for (size_t generator = 0; generator < abstractions_by_generator.size(); ++generator) {
        int number_affected_abs = abstractions_by_generator[generator].size();
        if (number_affected_abs < smallest_generator_size) {
            smallest_generator_size = number_affected_abs;
            smallest_generator_index = generator;
        }
    }
    if (smallest_generator_index == -1) {
        cerr << "find smallest generator failed: no abstractions affected at all" << endl;
        exit_with(EXIT_CRITICAL_ERROR);
    }
    abs_to_merge.insert(abstractions_by_generator[smallest_generator_index].begin(),
                        abstractions_by_generator[smallest_generator_index].end());
}

bool Symmetries::find_symmetries(const vector<Abstraction *>& abstractions,
                                 vector<set<int> > &affected_abstractions_by_generator,
                                 vector<set<int> > &mapped_abstractions_by_generator,
                                 vector<set<int> > &affected_not_mapped_abs_by_gen,
                                 vector<vector<vector<int> > > &cycles_by_generator,
                                 vector<int> &atomic_generators,
                                 vector<int> &local_generators) {
    // Find (non) abstraction stabilized symmetries for abstractions depending
    // on the chosen option.
    // When returning, affected_abstractions_by_generator contains the
    // abstraction indices that are affected by each generator,
    // atomic_generators contains the indices of atomic generators, and
    // local_generators contains the indices of local generators.
    // Returns true if any symmetry is found at all and false otherwise.

    // affected not mapped is the set difference of affected abs minus
    // mapped abs. There may exist generators both mapping abstractions onto
    // each other and affecting each of them individually as well.
    assert(affected_abstractions_by_generator.empty());
    assert(mapped_abstractions_by_generator.empty());
    assert(affected_not_mapped_abs_by_gen.empty());
    assert(cycles_by_generator.empty());

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
    affected_not_mapped_abs_by_gen.resize(num_generators, set<int>());
    cycles_by_generator.resize(num_generators, vector<vector<int> >());
    for (int gen_index = 0; gen_index < num_generators; ++gen_index) {
        cout << "generator " << gen_index << endl;
        set<int> &affected_abs = affected_abstractions_by_generator[gen_index];
        set<int> &mapped_abs = mapped_abstractions_by_generator[gen_index];
        set<int> &affected_not_mapped = affected_not_mapped_abs_by_gen[gen_index];
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
            vector<vector<int> > &cycles = cycles_by_generator[gen_index];
            vector<bool> marked(num_abstractions, false);
            for (size_t abs_index = 0; abs_index < num_abstractions; ++abs_index) {
                if (!marked[abs_index] && abstractions[abs_index]) {
                    marked[abs_index] = true;
                    unsigned int to_index = get_generator(gen_index)->get_value(abs_index);
                    if (to_index != abs_index) {
                        int from_index = abs_index;
                        vector<int> cycle;
                        cycle.push_back(from_index);
                        while (to_index != abs_index) {
                            marked[to_index] = true;
                            cycle.push_back(to_index);
                            from_index = to_index;
                            to_index = get_generator(gen_index)->get_value(from_index);
                        }
                        cycles.push_back(cycle);
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
        set_difference(affected_abs.begin(), affected_abs.end(),
                       mapped_abs.begin(), mapped_abs.end(),
                       inserter(affected_not_mapped, affected_not_mapped.begin()));
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
        cout << "Abstracting " << abstractions[i]->description() << endl;
        abstractions[i]->dump();
        abstractions[i]->dump_state();
        const EquivRel &equiv_rel = equivalence_relations[i];
        for (size_t j = 0; j < equiv_rel.size(); ++j) {
            cout << "class " << j << endl;
            const EquivClass &equiv_class = equiv_rel[j];
            for (EquivClass::const_iterator it = equiv_class.begin();
                 it != equiv_class.end(); ++it) {
                cout << *it << " ";
            }
            cout << endl;
        }
        abstractions[i]->apply_abstraction(equivalence_relations[i]);
    }
    cout << "Done abstracting. [t=" << g_timer << "]" << endl;
    cout << "==========================================================================================" << endl;
}

const Permutation* Symmetries::get_generator(int ind) const {
    assert(ind >= 0 && ind < get_num_generators());
    return gc.get_generators()[ind];
}
