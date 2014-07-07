#include "symmetries.h"

#include "scc.h"
#include "symmetry_generator.h"

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
    find_symmetries(abstractions);
    int smallest_generator_affected_abstractions_index = -1;
    int smallest_generator_affected_abstrations_size = numeric_limits<int>::max();
    int smallest_local_generator_mapped_abstractions_index = -1;
    int smallest_local_generator_mapped_abstractions__size = numeric_limits<int>::max();
    //int largest_atomic_cycle_index = -1;
    //int largest_atomic_cycle_size = 0;
    //int largest_local_cycle_index = -1;
    //int largest_local_cycle_size = 0;
    //vector<int> atomic_cycles;
    //vector<int> local_cycles;
    vector<int> atomic_generators;
    vector<int> local_generators;

    for (size_t generator_index = 0; generator_index < get_num_generators(); ++generator_index) {
        const SymmetryGenerator *generator = get_symmetry_generator(generator_index);
        const vector<int> &internally_affected_abstractions = generator->get_internally_affected_abstractions();
        const vector<int> &mapped_abstractions = generator->get_mapped_abstractions();
        const vector<int> &overall_affected_abstractions = generator->get_overall_affected_abstractions();
        //const vector<vector<int> > &cycles = generator->get_cycles();

        int number_overall_affected_abstractions = overall_affected_abstractions.size();
        if (number_overall_affected_abstractions < 1) {
            cerr << "Something is wrong! The generator is the identity generator." << endl;
            exit_with(EXIT_CRITICAL_ERROR);
        }
        if (number_overall_affected_abstractions == 1) {
            atomic_generators.push_back(generator_index);
        } else {
            if (number_overall_affected_abstractions < smallest_generator_affected_abstrations_size) {
                smallest_generator_affected_abstrations_size = number_overall_affected_abstractions;
                smallest_generator_affected_abstractions_index = generator_index;
            }
        }

        int number_mapped_abstractions = mapped_abstractions.size();
        if (number_mapped_abstractions == 0) {
            local_generators.push_back(generator_index);
        } else {
            if (number_mapped_abstractions < smallest_local_generator_mapped_abstractions__size) {
                smallest_local_generator_mapped_abstractions__size = number_mapped_abstractions;
                smallest_local_generator_mapped_abstractions_index = generator_index;
            }
        }

//            if (cycles.size() > 0) {
//                local_cycles.push_back(generator_index);
//                /*int cycles_size = 0;
//                for (size_t i = 0; i < cycles.size(); ++i) {
//                    cycles_size += cycles[i].size();
//                }
//                if (cycles_size > largest_local_cycle_size) {
//                    largest_local_cycle_size = cycles_size;
//                    largest_local_cycle_index = generator_index;
//                }*/
//                if (cycles.size() == 1) {
//                    vector<int> affected_minus_mapped;
//                    // affected minues mapped is the set difference of affected abs minus
//                    // mapped abs. If it is non-empty, there are abstractions not mapped
//                    // which are affected by the generator. As a result, the generator
//                    // is not atomic.
//                    // NOTE: No generator may both internally affect an abstraction
//                    // and map it onto another one at the same time, due to the
//                    // definition of the PDG.
//                    set_difference(internally_affected_abstractions.begin(), internally_affected_abstractions.end(),
//                                   mapped_abstractions.begin(), mapped_abstractions.end(),
//                                   inserter(affected_minus_mapped, affected_minus_mapped.begin()));

//                    if (affected_minus_mapped.empty()) {
//                        atomic_cycles.push_back(generator_index);
//                        //if (cycles_size > largest_atomic_cycle_size) {
//                        //    largest_atomic_cycle_size = cycles_size;
//                        //    largest_atomic_cycle_index = generator_index;
//                        //}
//                    }
//                }
//            }

        switch (number_overall_affected_abstractions) {
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

        // dump generator properties
        cout << "Generator " << generator_index << endl;
        for (size_t i = 0; i < mapped_abstractions.size(); ++i) {
            int abs_index = mapped_abstractions[i];
            int to_index = generator->get_value(abs_index);
            cout << abstractions[abs_index]->description() << " mapped to " <<
                    abstractions[to_index]->description();
            if (generator->internally_affects(abs_index))
                cout << " (and also internally affected)";
            cout << endl;
        }
        for (size_t i = 0; i < internally_affected_abstractions.size(); ++i) {
            int abs_index = internally_affected_abstractions[i];
            assert(!generator->maps(abs_index));
            cout << abstractions[abs_index]->description() << " internally affected" << endl;
        }
    }

    switch (type_of_symmetries) {
        case ATOMIC: {
            if (!atomic_generators.empty()) {
                // apply atomic symmetries
                ++number_of_applied_symmetries;
                cout << "Applying all atomic symmetries" << endl;
                apply_symmetries(abstractions, atomic_generators);
            }

            // NOTE: atomic generators and atomic cycles are always disjoint
            // sets of generators, as an atomic generator affects exactly
            // one abstraction, whereas an atomic cycle consists of at least
            // two abstractions.
            //if (atomic_generators.empty()/* && atomic_cycles.empty()*/) {
            if (smallest_generator_affected_abstractions_index != -1) {
                // use the "smallest" symmetry to merge for
                const vector<int> &overall_affected_abstractions =
                        get_symmetry_generator(smallest_generator_affected_abstractions_index)->
                        get_overall_affected_abstractions();
                abs_to_merge.insert(overall_affected_abstractions.begin(), overall_affected_abstractions.end());
                //return make_pair(number_of_applied_symmetries, number_of_collapsed_abstractions);
            }

//                // copy all cycles into a new data structure
//                vector<vector<int> > new_cycles;
//                for (size_t i = 0; i < atomic_cycles.size(); ++i) {
//                    int generator_index = atomic_cycles[i];
//                    const SymmetryGenerator *generator = get_symmetry_generator(generator_index);
//                    const vector<vector<int> > &cycles = generator->get_cycles();
//                    assert(cycles.size() == 1);
//                    for (size_t cycle_no = 0; cycle_no < cycles.size(); ++cycle_no) {
//                        const vector<int> &collapsed_abs = cycles[cycle_no];
//                        new_cycles.push_back(vector<int>(collapsed_abs));
//                    }
//                }

//                // find and combine all "overlapping" cycles
//                for (int i = 0; i < new_cycles.size() - 1; ++i) {
//                    vector<int> &cycle1 = new_cycles[i];
//                    //cout << "considering cycle " << cycle1 << endl;
//                    for (int j = i + 1; j < new_cycles.size(); ++j) {
//                        const vector<int> &cycle2 = new_cycles[j];
//                        //cout << "and " << cycle2 << endl;

//                        vector<int> intersection;
//                        set_intersection(cycle1.begin(), cycle1.end(),
//                                         cycle2.begin(), cycle2.end(),
//                                         back_inserter(intersection));
//                        if (!intersection.empty()) {
//                            // cycles overlap
//                            vector<int> unified_cycle;
//                            set_union(cycle1.begin(), cycle1.end(),
//                                      cycle2.begin(), cycle2.end(),
//                                      back_inserter(unified_cycle));
//                            // replace cycle1 by the union of cycle1 and cycle2
//                            cycle1.swap(unified_cycle);
//                            // erase cycle2
//                            new_cycles.erase(new_cycles.begin() + j);
//                            // repeat outer loop
//                            --i;

//                            //cout << "overlap!" << endl;
//                            //for (size_t x = 0; x < new_cycles.size(); ++x)
//                            //    cout << new_cycles[x] << endl;
//                            break;
//                        }
//                    }
//                }

//                // apply all cycles
//                vector<bool> mapped_abstractions(abstractions.size(), false);
//                abstractions[7]->dump();
//                for (size_t cycle_no = 0; cycle_no < new_cycles.size(); ++cycle_no) {
//                    const vector<int> &cycle = new_cycles[cycle_no];
//                    cout << "Collapsing cycle: " << cycle << endl;
//                    Abstraction *chosen_representative = abstractions[cycle[0]];
//                    if (!mapped_abstractions[cycle[0]]) {
//                        mapped_abstractions[cycle[0]] = true;
//                    } else {
//                        cerr << "Abstraction is included in several cycles" << endl;
//                        exit_with(EXIT_CRITICAL_ERROR);
//                    }
//                    for (size_t j = 1; j < cycle.size(); ++j) {
//                        size_t abs_index = cycle[j];
//                        if (!mapped_abstractions[abs_index]) {
//                            mapped_abstractions[abs_index] = true;
//                        } else {
//                            cerr << "Abstraction is included in several cycles" << endl;
//                            exit_with(EXIT_CRITICAL_ERROR);
//                        }
//                        Abstraction *abs = abstractions[abs_index];
//                        chosen_representative->merge_abstraction_into(abs);
//                        chosen_representative->normalize();
//                        abs->release_memory();
//                        abstractions[abs_index] = 0;
//                    }
//                    number_of_collapsed_abstractions += cycle.size() - 1;
//                }
            break;
        }
        case LOCAL: {
            if (!local_generators.empty()) {
                // apply local symmetries
                ++number_of_applied_symmetries;
                cout << "Applying all local symmetries" << endl;
                apply_symmetries(abstractions, local_generators);
            }

            //if (local_generators.empty()) {
            if (smallest_local_generator_mapped_abstractions_index != -1) {
                // use the "smallest" symmetry to merge for
                const vector<int> &mapped_abstractions =
                        get_symmetry_generator(smallest_local_generator_mapped_abstractions_index)->
                        get_mapped_abstractions();
                abs_to_merge.insert(mapped_abstractions.begin(), mapped_abstractions.end());
                //return make_pair(number_of_applied_symmetries, number_of_collapsed_abstractions);
            }
            break;
        }
        case MERGE_ONLY: {
            if (smallest_generator_affected_abstractions_index != -1) {
                const vector<int> &overall_affected_abstractions =
                        get_symmetry_generator(smallest_generator_affected_abstractions_index)->
                        get_overall_affected_abstractions();
                abs_to_merge.insert(overall_affected_abstractions.begin(), overall_affected_abstractions.end());
            }
        }
    }
    return make_pair(number_of_applied_symmetries, number_of_collapsed_abstractions);
}

bool Symmetries::find_symmetries(const vector<Abstraction *>& abstractions) {
    // Find (non) abstraction stabilized symmetries for abstractions depending
    // on the chosen option.
    // Returns true if any symmetry is found at all and false otherwise.

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
    return true;
}

void Symmetries::apply_symmetries(const vector<Abstraction *> &abstractions,
                                  const vector<int> &generator_indices) const {
    if (get_num_generators() == 0) {
        cerr << "You first have to find symmetries before you can apply one of them!" << endl;
        exit_with(EXIT_CRITICAL_ERROR);
    }
    cout << "Creating equivalence relations from symmetries. [t=" << g_timer << "]" << endl;

    // TODO: this does not work for abstractions which are entirely mapped to
    // another abstractions (see cycles above)! We need to make sure here that
    // states of an abstraction are only mapped to states of the same abstraction

    // Abstracting by the generators as follows:
    // Creating a graph with nodes being the abstract states and the edges represent the connections as given by the generators.
    // It seems like this process should better be done in all abstractions in parallel,
    // since we could exploit all the compositions of these abstractions this way as well.
    // Later we can reduce the first part - the abstraction vars and save some space/time
    unsigned int num_states = get_sym_gen_info().num_abs_and_states;
    unsigned int num_abstractions = get_sym_gen_info().num_abstractions;
    // The graph is represented by vector of to_nodes for each node. (Change to sets?)
    vector<vector<unsigned int> > graph(num_states, vector<unsigned int>());
    for (unsigned int index = num_abstractions; index < get_sym_gen_info().num_abs_and_states; ++index) {
        int abs_index = get_sym_gen_info().get_var_by_index(index);
        if (abstractions[abs_index]) {
            for (size_t i = 0; i < generator_indices.size(); ++i) {
                // Going over the generators, for each just add the edges.
                unsigned int to_index = get_symmetry_generator(generator_indices[i])->get_value(index);
                if (index != to_index)
                    graph[index].push_back(to_index);
            }
        }
    }
    SCC scc(graph);
    const vector<vector<unsigned int> >& result = scc.get_result();

    // Generate final result. Going over the result, putting the nodes to their respective places.
    vector<EquivRel> equivalence_relations(abstractions.size(), EquivRel());
    for (size_t abs_index = 0; abs_index < abstractions.size(); ++abs_index) {
        if (abstractions[abs_index])
            equivalence_relations[abs_index].resize(result.size());
    }
    for (unsigned int eqiv=0; eqiv < result.size(); eqiv++) {
        for (unsigned int i=0; i < result[eqiv].size(); i++) {
            unsigned int idx = result[eqiv][i];
            if (idx < get_sym_gen_info().num_abstractions)
                continue;
            pair<int, AbstractStateRef> vals = get_sym_gen_info().get_var_val_by_index(idx);
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
        /*abstractions[i]->dump();
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
        }*/
        abstractions[i]->apply_abstraction(equivalence_relations[i]);
    }
    cout << "Done abstracting. [t=" << g_timer << "]" << endl;
    cout << "==========================================================================================" << endl;
}

const SymmetryGenerator* Symmetries::get_symmetry_generator(int ind) const {
    assert(ind >= 0 && ind < get_num_generators());
    return gc.get_symmetry_generators()[ind];
}
