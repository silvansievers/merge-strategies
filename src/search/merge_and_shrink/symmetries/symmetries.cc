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
      symmetries_for_shrinking(SymmetriesForShrinking(options.get_enum("symmetries_for_shrinking"))),
      symmetries_for_merging(SymmetriesForMerging(options.get_enum("symmetries_for_merging"))),
      internal_merging(InternalMerging(options.get_enum("internal_merging"))),
      build_stabilized_pdg(options.get<bool>("build_stabilized_pdg")) {
}

pair<int, int> Symmetries::find_and_apply_symmetries(vector<Abstraction *> &abstractions,
                                                     vector<vector<int> > &abs_to_merge) {
    assert(abs_to_merge.empty());
    int number_of_applied_symmetries = 0;
    int number_of_collapsed_abstractions = 0;
    find_symmetries(abstractions);
    int chosen_generator_affected_abstractions_index = -1;
    int chosen_generator_mapped_abstractions_index = -1;
    int smallest_generator_affected_abstrations_size = numeric_limits<int>::max();
    int smallest_generator_mapped_abstractions_size = numeric_limits<int>::max();
    int largest_generator_affected_abstrations_size = 0;
    int largest_generator_mapped_abstractions_size = 0;
    vector<int> atomic_generators;
    vector<int> local_generators;

    for (size_t generator_index = 0; generator_index < get_num_generators(); ++generator_index) {
        const SymmetryGenerator *generator = get_symmetry_generator(generator_index);
        const vector<int> &internally_affected_abstractions = generator->get_internally_affected_abstractions();
        const vector<int> &mapped_abstractions = generator->get_mapped_abstractions();
        const vector<int> &overall_affected_abstractions = generator->get_overall_affected_abstractions();

        int number_overall_affected_abstractions = overall_affected_abstractions.size();
        if (number_overall_affected_abstractions < 1) {
            cerr << "Something is wrong! The generator is the identity generator." << endl;
            exit_with(EXIT_CRITICAL_ERROR);
        }
        if (number_overall_affected_abstractions == 1) {
            atomic_generators.push_back(generator_index);
        } else {
            if (symmetries_for_merging == SMALLEST
                    && number_overall_affected_abstractions < smallest_generator_affected_abstrations_size) {
                smallest_generator_affected_abstrations_size = number_overall_affected_abstractions;
                chosen_generator_affected_abstractions_index = generator_index;
            }
            if (symmetries_for_merging == LARGEST
                    && number_overall_affected_abstractions > largest_generator_affected_abstrations_size) {
                largest_generator_affected_abstrations_size = number_overall_affected_abstractions;
                chosen_generator_affected_abstractions_index = generator_index;
            }
        }

        int number_mapped_abstractions = mapped_abstractions.size();
        if (number_mapped_abstractions == 0) {
            local_generators.push_back(generator_index);
        } else {
            if (symmetries_for_merging == SMALLEST
                    && number_mapped_abstractions < smallest_generator_mapped_abstractions_size) {
                smallest_generator_mapped_abstractions_size = number_mapped_abstractions;
                chosen_generator_mapped_abstractions_index = generator_index;
            }
            if (symmetries_for_merging == LARGEST
                    && number_mapped_abstractions > largest_generator_mapped_abstractions_size) {
                largest_generator_mapped_abstractions_size = number_mapped_abstractions;
                chosen_generator_mapped_abstractions_index = generator_index;
            }
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

    switch (symmetries_for_shrinking) {
        case ATOMIC: {
            if (!atomic_generators.empty()) {
                // apply atomic symmetries
                ++number_of_applied_symmetries;
                cout << "Applying all atomic symmetries" << endl;
                apply_symmetries(abstractions, atomic_generators);
            }
            break;
        }
        case LOCAL: {
            if (!local_generators.empty()) {
                // apply local symmetries
                ++number_of_applied_symmetries;
                cout << "Applying all local symmetries" << endl;
                apply_symmetries(abstractions, local_generators);
            }
            break;
        }
        case NONE: {
            break;
        }
    }

    if (symmetries_for_merging != ZERO) {
        if (symmetries_for_shrinking == ATOMIC || symmetries_for_shrinking == NONE) {
            if (chosen_generator_affected_abstractions_index != -1) {
                const vector<int> &overall_affected_abstractions =
                        get_symmetry_generator(chosen_generator_affected_abstractions_index)->
                        get_overall_affected_abstractions();
                abs_to_merge.push_back(vector<int>(overall_affected_abstractions));
            }
        } else if (symmetries_for_shrinking == LOCAL) {
            if (chosen_generator_mapped_abstractions_index != -1) {
                if (internal_merging == LINEAR) {
                    const vector<int> &mapped_abstractions =
                            get_symmetry_generator(chosen_generator_mapped_abstractions_index)->
                            get_mapped_abstractions();
                    abs_to_merge.push_back(vector<int>(mapped_abstractions));
                } else if (internal_merging == NON_LINEAR) {
                    get_symmetry_generator(chosen_generator_mapped_abstractions_index)->compute_cycles(abs_to_merge);
                }
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
                if (get_symmetry_generator(generator_indices[i])->get_value(abs_index) == abs_index) {
                    // we only add an edge if the corresponding states belong to
                    // the same abstractions. in other words, we do not compute
                    // equivalence relations for mappings of abstractions, as
                    // these are not applied anyways.
                    unsigned int to_index = get_symmetry_generator(generator_indices[i])->get_value(index);
                    if (index != to_index)
                        graph[index].push_back(to_index);
                }
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
