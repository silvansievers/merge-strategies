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
    while (true) {
        //for (size_t i = 0; i < abstractions.size(); ++i) {
        //    cout << i << " " << abstractions[i] << endl;
        //}

        bool found_symmetry = find_symmetries(abstractions);
        if (found_symmetry) {
            int smallest_generator_affected_abstractions_index = -1;
            int smallest_generator_affected_abstrations_size = numeric_limits<int>::max();
            //int smallest_generator_mapped_abstractions_index = -1;
            //int smallest_generator_mapped_abstrations_size = numeric_limits<int>::max();
            int largest_atomic_cycle_index = -1;
            int largest_atomic_cycle_size = 0;
            int largest_local_cycle_index = -1;
            int largest_local_cycle_size = 0;
            vector<int> atomic_generators;
            vector<int> local_generators;

            for (size_t generator_index = 0; generator_index < get_num_generators(); ++generator_index) {
                const SymmetryGenerator *generator = get_symmetry_generator(generator_index);
                const vector<int> &internally_affected_abstractions = generator->get_internally_affected_abstractions();
                const vector<int> &mapped_abstractions = generator->get_mapped_abstractions();
                const vector<int> &overall_affected_abstractions = generator->get_overall_affected_abstractions();
                const vector<vector<int> > &cycles = generator->get_cycles();

                //int number_mapped_abs = mapped_abstractions.size();
                int number_overall_affected_abstractions = overall_affected_abstractions.size();

                if (number_overall_affected_abstractions < 1) {
                    cerr << "Something is wrong! The generator is the identity generator." << endl;
                    exit_with(EXIT_CRITICAL_ERROR);
                }
                if (number_overall_affected_abstractions == 1) {
                    atomic_generators.push_back(generator_index);
                    // every atomic symmetry is also a local symmetry
                    //assert(number_mapped_abs == 0);
                }

                //if (number_mapped_abs == 0) {
                    // local currently means no cyles, but neither atomic
                    local_generators.push_back(generator_index);
                //}

                if (number_overall_affected_abstractions < smallest_generator_affected_abstrations_size) {
                    smallest_generator_affected_abstrations_size = number_overall_affected_abstractions;
                    smallest_generator_affected_abstractions_index = generator_index;
                }
                /*if (number_mapped_abs < smallest_generator_mapped_abstrations_size) {
                    smallest_generator_mapped_abstrations_size = number_mapped_abs;
                    smallest_generator_mapped_abstractions_index = generator_index;
                }*/

                vector<int> intersection_internally_affected_mapped;
                // intersection affected mapped is the set intersection of affected abs and
                // mapped abs. There may exist generators both mapping abstractions onto
                // each other and affecting each of them individually as well, which
                // means that the generator is not atomic.
                // TODO: I think this can never happen
                set_intersection(internally_affected_abstractions.begin(), internally_affected_abstractions.end(),
                                 mapped_abstractions.begin(), mapped_abstractions.end(),
                                 inserter(intersection_internally_affected_mapped, intersection_internally_affected_mapped.begin()));
                if (!intersection_internally_affected_mapped.empty()) {
                    cerr << "Abstraction is mapped and internally affected!" << endl;
                    exit_with(EXIT_CRITICAL_ERROR);
                }

                if (cycles.size() > 0) {
                    int cycles_size = 0;
                    for (size_t i = 0; i < cycles.size(); ++i) {
                        cycles_size += cycles[i].size();
                    }
                    if (cycles_size > largest_local_cycle_size) {
                        largest_local_cycle_size = cycles_size;
                        largest_local_cycle_index = generator_index;
                    }
                    if (cycles.size() == 1) {
                        vector<int> affected_minus_mapped;
                        // affected minues mapped is the set difference of affected abs minus
                        // mapped abs. If it is non-empty, there are abstractions not mapped
                        // which are affected by the generator.
                        set_difference(internally_affected_abstractions.begin(), internally_affected_abstractions.end(),
                                       mapped_abstractions.begin(), mapped_abstractions.end(),
                                       inserter(affected_minus_mapped, affected_minus_mapped.begin()));

                        if (affected_minus_mapped.empty() && intersection_internally_affected_mapped.empty()) {
                            if (cycles_size > largest_atomic_cycle_size) {
                                largest_atomic_cycle_size = cycles_size;
                                largest_atomic_cycle_index = generator_index;
                            }
                        }
                    }
                }

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
                    if (!generator->maps(abs_index))
                        cout << abstractions[abs_index]->description() << " internally affected" << endl;
                }
            }

            switch (type_of_symmetries) {
            case ATOMIC:
                if (largest_atomic_cycle_index != -1) {
                    // there is a generator with exactly one cyclic mapping
                    ++number_of_applied_symmetries;
                    cout << "'Merging' abstractions by removing all but one abstraction from a cycle" << endl;
                    const SymmetryGenerator *generator = get_symmetry_generator(largest_atomic_cycle_index);
                    const vector<vector<int> > &cycles = generator->get_cycles();
                    assert(cycles.size() == 1);
                    const vector<int> &collapsed_abs = cycles[0];
                    Abstraction *chosen_representative = abstractions[collapsed_abs[0]];
                    for (size_t i = 1; i < collapsed_abs.size(); ++i) {
                        size_t abs_index = collapsed_abs[i];
                        Abstraction *abs = abstractions[abs_index];
                        chosen_representative->merge_abstraction_into(abs);
                        chosen_representative->normalize();
                        abs->release_memory();
                        abstractions[abs_index] = 0;
                    }
                    number_of_collapsed_abstractions += collapsed_abs.size() - 1;
                } else {
                    // no atomic cyclic mappings
                    if (atomic_generators.empty()) {
                        // find the "smallest" symmetry to merge for
                        assert(smallest_generator_affected_abstractions_index != -1);
                        const vector<int> &overall_affected_abstractions =
                                get_symmetry_generator(smallest_generator_affected_abstractions_index)->
                                get_overall_affected_abstractions();
                        abs_to_merge.insert(overall_affected_abstractions.begin(), overall_affected_abstractions.end());
                        return make_pair(number_of_applied_symmetries, number_of_collapsed_abstractions);
                    } else {
                        // apply atomic symmetries
                        ++number_of_applied_symmetries;
                        apply_symmetries(abstractions, atomic_generators);
                    }
                }
                break;
            case LOCAL:
                if (largest_local_cycle_index != -1) {
                    // there is a generator with cyclic abstraction mappings

                    // first, apply the generator to the set of all abstractions.
                    // this is required as to make sure that the internally
                    // affected abstractions by the generator are correctly
                    // modified.
                    // TODO: only do this for generators that actually affect
                    // abstractions internally besides having cyclic mappings
                    // -> this information should be computed in the generators
                    // already.
                    vector<int> generator_to_be_applied;
                    generator_to_be_applied.push_back(largest_local_cycle_index);
                    apply_symmetries(abstractions, generator_to_be_applied);

                    // second, for all cycles, choose a single representing
                    // abstraction
                    ++number_of_applied_symmetries;
                    cout << "'Merging' abstractions by removing all but one abstraction from (a) cycle(s)" << endl;
                    const SymmetryGenerator *generator = get_symmetry_generator(largest_local_cycle_index);
                    const vector<vector<int> > &cycles = generator->get_cycles();
                    for (size_t cycle_no = 0; cycle_no < cycles.size(); ++cycle_no) {
                        const vector<int> &collapsed_abs = cycles[cycle_no];
                        Abstraction *chosen_representative = abstractions[collapsed_abs[0]];
                        for (size_t i = 1; i < collapsed_abs.size(); ++i) {
                            size_t abs_index = collapsed_abs[i];
                            Abstraction *abs = abstractions[abs_index];
                            chosen_representative->merge_abstraction_into(abs);
                            chosen_representative->normalize();
                            abs->release_memory();
                            abstractions[abs_index] = 0;
                        }
                        number_of_collapsed_abstractions += collapsed_abs.size() - 1;
                    }
                } else {
                    // no cyclic mappings at all
                    ++number_of_applied_symmetries;
                    apply_symmetries(abstractions, local_generators);
                    /*if (local_generators.empty() && atomic_generators.empty()) {
                        if (smallest_generator_mapped_abstractions_index == -1) {
                            cerr << "No atomic or local generators found, but mapped abstractions "
                                    "is also empty!" << endl;
                            exit_with(EXIT_CRITICAL_ERROR);
                        }
                        const vector<int> &mapped_abstractions =
                                get_symmetry_generator(smallest_generator_mapped_abstractions_index)->
                                get_mapped_abstractions();
                        abs_to_merge.insert(mapped_abstractions.begin(), mapped_abstractions.end());
                        return make_pair(number_of_applied_symmetries, 0);
                    } else {

                    }*/
                    break;
                }
            }
        } else {
            break;
        }
    }
    return make_pair(number_of_applied_symmetries, number_of_collapsed_abstractions);
}

bool Symmetries::find_symmetries(const vector<Abstraction *>& abstractions) {
    // Find (non) abstraction stabilized symmetries for abstractions depending
    // on the chosen option.
    // When returning, affected_abstractions_by_generator contains the
    // abstraction indices that are affected by each generator,
    // atomic_generators contains the indices of atomic generators, and
    // local_generators contains the indices of local generators.
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
