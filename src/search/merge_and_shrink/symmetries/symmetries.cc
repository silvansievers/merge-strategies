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

bool Symmetries::find_and_apply_symmetries(vector<Abstraction *> &abstractions,
                                           vector<pair<int, int> > &merge_order) {
    find_symmetries(abstractions);
    int chosen_generator_for_merging = -1;
    int smallest_generator_affected_abstrations_size = numeric_limits<int>::max();
    int smallest_generator_mapped_abstractions_size = numeric_limits<int>::max();
    int largest_generator_affected_abstrations_size = 0;
    int largest_generator_mapped_abstractions_size = 0;
    vector<int> atomic_generators;
    vector<int> local_generators;

    // go over the generators and classify them into atomic, local or general
    // ones. also store information about the smallest/largest generator
    // with respect to overall affected abstractions or mapped abstractions,
    // depending on the chosen setting.
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
            if (symmetries_for_merging == LEAST_OVERALL_AFFECTED
                    && number_overall_affected_abstractions < smallest_generator_affected_abstrations_size) {
                smallest_generator_affected_abstrations_size = number_overall_affected_abstractions;
                chosen_generator_for_merging = generator_index;
            }
            if (symmetries_for_merging == MOST_OVERALL_AFFECTED
                    && number_overall_affected_abstractions > largest_generator_affected_abstrations_size) {
                largest_generator_affected_abstrations_size = number_overall_affected_abstractions;
                chosen_generator_for_merging = generator_index;
            }
        }

        int number_mapped_abstractions = mapped_abstractions.size();
        if (number_mapped_abstractions == 0) {
            local_generators.push_back(generator_index);
        } else {
            if (symmetries_for_merging == LEAST_MAPPED
                    && number_mapped_abstractions < smallest_generator_mapped_abstractions_size) {
                smallest_generator_mapped_abstractions_size = number_mapped_abstractions;
                chosen_generator_for_merging = generator_index;
            }
            if (symmetries_for_merging == MOST_MAPPED
                    && number_mapped_abstractions > largest_generator_mapped_abstractions_size) {
                largest_generator_mapped_abstractions_size = number_mapped_abstractions;
                chosen_generator_for_merging = generator_index;
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

    // apply symmetries if possible
    bool applied_symmetries = false;
    if (symmetries_for_shrinking == ATOMIC && !atomic_generators.empty()) {
        // apply atomic symmetries
        cout << "Applying all atomic symmetries" << endl;
        apply_symmetries(abstractions, atomic_generators);
        applied_symmetries = true;
    } else if (symmetries_for_shrinking == LOCAL && !local_generators.empty()) {
        // apply local symmetries
        cout << "Applying all local symmetries" << endl;
        apply_symmetries(abstractions, local_generators);
        applied_symmetries = true;
    }

    if (symmetries_for_merging != NO_MERGING && chosen_generator_for_merging != -1) {
        vector<vector<int> > merge_non_linear_abstractions;
        vector<int> merge_linear_abstractions;
        const SymmetryGenerator *generator =
                get_symmetry_generator(chosen_generator_for_merging);

        // Always include all mapped abstractions
        if (internal_merging == NON_LINEAR
                || internal_merging == NON_LINEAR_INCOMPLETE) {
            // if the internal merge strategy is non linear, we need to
            // compute the actual cycles of abstraction mappings.
            generator->compute_cycles(merge_non_linear_abstractions);
        } else if (internal_merging == LINEAR) {
            // if the internal merge strategy is linear, we simply collect
            // all mapped abstractions (i.e. the same abstractions as above,
            // but we do not compute cycle information)
            const vector<int> &mapped_abstractions =
                    generator->get_mapped_abstractions();
            merge_linear_abstractions.insert(merge_linear_abstractions.end(),
                                             mapped_abstractions.begin(),
                                             mapped_abstractions.end());
        }

        // If merging for least/most number of overall affected abstactions,
        // also include the non-mapped, i.e. internally affected abstractions
        // (always as to be linearly merged abstractions)
        if (symmetries_for_merging == LEAST_OVERALL_AFFECTED
                || symmetries_for_merging == MOST_OVERALL_AFFECTED) {
            const vector<int> &internally_affected_abstractions =
                    generator->get_internally_affected_abstractions();
            merge_linear_abstractions.insert(merge_linear_abstractions.end(),
                                             internally_affected_abstractions.begin(),
                                             internally_affected_abstractions.end());
        }

        // compute a merge tree
        assert(merge_order.empty());
        int number_of_abstractions = abstractions.size();
        int number_of_merges = 0;
        vector<int> merge_linear_indices;
        for (size_t cycle_no = 0; cycle_no < merge_non_linear_abstractions.size(); ++cycle_no) {
            // go over the cycles and compute a non-linear merge order.
            const vector<int> &cycle = merge_non_linear_abstractions[cycle_no];
            size_t abs_index_1 = cycle[0];
            for (size_t i = 1; i < cycle.size(); ++i) {
                size_t abs_index_2 = cycle[i];
                merge_order.push_back(make_pair(abs_index_1, abs_index_2));
                abs_index_1 = number_of_abstractions + number_of_merges;
                ++number_of_merges;
            }
            // number_of_abstractions + number_of_merges always is the *next*
            // position where a new merged abstraction will be stored at.
            // here, we need the *last* position where the abstraction
            // resulting from merging the cycle was stored, hence the -1.
            merge_linear_indices.push_back(number_of_abstractions + number_of_merges - 1);
        }

        if (internal_merging != NON_LINEAR_INCOMPLETE) {
            // merge_linear_indices possibly contains abstractions that have been
            // non-linearly merged from information about cycles.
            // here we add abstractions that need to be merged linearly anyways
            merge_linear_indices.insert(merge_linear_indices.end(),
                                        merge_linear_abstractions.begin(),
                                        merge_linear_abstractions.end());

            // go over all abstractions that (now) need to be merged linearly
            size_t abs_index_1 = merge_linear_indices[0];
            for (size_t i = 1; i < merge_linear_indices.size(); ++i) {
                size_t abs_index_2 = merge_linear_indices[i];
                merge_order.push_back(make_pair(abs_index_1, abs_index_2));
                abs_index_1 = number_of_abstractions + number_of_merges;
                ++number_of_merges;
            }
        }

        cout << "current number of abstractions " << number_of_abstractions << endl;
        cout << "chosen internal merge order: " << endl;
        for (size_t i = 0; i < merge_order.size(); ++i) {
            cout << merge_order[i].first << ", " << merge_order[i].second << endl;
        }
    }

    return applied_symmetries;
}

void Symmetries::find_symmetries(const vector<Abstraction *>& abstractions) {
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
    cout << "Done searching for symmetries. [t=" << g_timer << "]" << endl;
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
