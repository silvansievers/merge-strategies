#include "symmetries.h"

#include "scc.h"
#include "symmetry_generator.h"

#include "../transition_system.h"

#include "../../globals.h"
#include "../../option_parser.h"
#include "../../timer.h"
#include "../../utilities.h"

#include <cassert>
#include <iostream>
#include <limits>

using namespace std;

Symmetries::Symmetries(const Options &options)
    : gc(options),
      symmetries_for_shrinking(SymmetriesForShrinking(options.get_enum("symmetries_for_shrinking"))),
      symmetries_for_merging(SymmetriesForMerging(options.get_enum("symmetries_for_merging"))),
      external_merging(ExternalMerging(options.get_enum("external_merging"))),
      internal_merging(InternalMerging(options.get_enum("internal_merging"))),
      bliss_time(0) {
}

bool Symmetries::find_and_apply_symmetries(const vector<TransitionSystem *> &transition_systems,
                                           vector<pair<int, int> > &merge_order) {
    bliss_time = gc.compute_generators(transition_systems);
    if (get_num_generators() == 0 || is_bliss_limit_reached()) {
        return false;
    }


    int chosen_generator_for_merging = -1;
    int smallest_generator_affected_transition_systems_size = numeric_limits<int>::max();
    int smallest_generator_mapped_transition_systems_size = numeric_limits<int>::max();
    int largest_generator_affected_transition_systems_size = 0;
    int largest_generator_mapped_transition_systems_size = 0;
    vector<int> atomic_generators;
    vector<int> local_generators;

    // go over the generators and classify them into atomic, local or general
    // ones. also store information about the smallest/largest generator
    // with respect to overall affected transition systems or mapped transition systems,
    // depending on the chosen setting.
    for (int generator_index = 0; generator_index < get_num_generators(); ++generator_index) {
        const SymmetryGenerator *generator = get_symmetry_generator(generator_index);
        const vector<int> &internally_affected_transition_systems = generator->get_internally_affected_transition_systems();
        const vector<int> &mapped_transition_systems = generator->get_mapped_transition_systems();
        const vector<int> &overall_affected_transition_systems = generator->get_overall_affected_transition_systems();

        int number_overall_affected_transition_systems = overall_affected_transition_systems.size();
        if (number_overall_affected_transition_systems < 1) {
            cerr << "Something is wrong! The generator is the identity generator." << endl;
            exit_with(EXIT_CRITICAL_ERROR);
        }
        if (number_overall_affected_transition_systems == 1) {
            atomic_generators.push_back(generator_index);
        } else {
            if (external_merging == MERGE_FOR_ATOMIC) {
                if (symmetries_for_merging == SMALLEST
                        && number_overall_affected_transition_systems
                        < smallest_generator_affected_transition_systems_size) {
                    smallest_generator_affected_transition_systems_size
                            = number_overall_affected_transition_systems;
                    chosen_generator_for_merging = generator_index;
                } else if (symmetries_for_merging == LARGEST
                              && number_overall_affected_transition_systems
                              > largest_generator_affected_transition_systems_size) {
                    largest_generator_affected_transition_systems_size
                            = number_overall_affected_transition_systems;
                    chosen_generator_for_merging = generator_index;
                }
            }
        }

        int number_mapped_transition_systems = mapped_transition_systems.size();
        if (number_mapped_transition_systems == 0) {
            // note that this also includes atomic generators
            local_generators.push_back(generator_index);
        } else {
            if (external_merging == MERGE_FOR_LOCAL) {
                if (symmetries_for_merging == SMALLEST
                        && number_mapped_transition_systems
                        < smallest_generator_mapped_transition_systems_size) {
                    smallest_generator_mapped_transition_systems_size
                            = number_mapped_transition_systems;
                    chosen_generator_for_merging = generator_index;
                }
                if (symmetries_for_merging == LARGEST
                        && number_mapped_transition_systems
                        > largest_generator_mapped_transition_systems_size) {
                    largest_generator_mapped_transition_systems_size
                            = number_mapped_transition_systems;
                    chosen_generator_for_merging = generator_index;
                }
            }
        }

        // dump generator properties
        cout << "Generator " << generator_index << endl;
        for (size_t i = 0; i < mapped_transition_systems.size(); ++i) {
            int abs_index = mapped_transition_systems[i];
            int to_index = generator->get_value(abs_index);
            cout << transition_systems[abs_index]->tag() << " mapped to " <<
                    transition_systems[to_index]->tag();
            if (generator->internally_affects(abs_index))
                cout << " (and also internally affected)";
            cout << endl;
        }
        for (size_t i = 0; i < internally_affected_transition_systems.size(); ++i) {
            int abs_index = internally_affected_transition_systems[i];
            assert(!generator->maps(abs_index));
            cout << transition_systems[abs_index]->tag() << " internally affected" << endl;
        }
    }

    // apply symmetries if possible
    bool applied_symmetries = false;
    if (symmetries_for_shrinking == ATOMIC && !atomic_generators.empty()) {
        // apply atomic symmetries
        cout << "Applying all atomic symmetries" << endl;
        apply_symmetries(transition_systems, atomic_generators);
        applied_symmetries = true;
    } else if ((symmetries_for_shrinking == LOCAL)
               && !local_generators.empty()) {
        // apply local symmetries
        cout << "Applying all local symmetries" << endl;
        apply_symmetries(transition_systems, local_generators);
        applied_symmetries = true;
    }

    if (symmetries_for_merging != NO_MERGING && chosen_generator_for_merging != -1) {
        vector<vector<int> > cycles;
        vector<int> merge_linear_transition_systems;
        const SymmetryGenerator *generator =
                get_symmetry_generator(chosen_generator_for_merging);

        // Always include all mapped transition systems
        if (internal_merging == NON_LINEAR
                || external_merging == MERGE_FOR_LOCAL) {
            // if the internal merge strategy is non linear or we only want
            // to merge every cycle (non linearly), we need to
            // compute the actual cycles of transition system mappings.
            generator->compute_cycles(cycles);
        } else if (internal_merging == LINEAR) {
            // if the internal merge strategy is linear, we simply collect
            // all mapped transition systems (i.e. the same transition systems as above,
            // but we do not compute cycle information)
            const vector<int> &mapped_transition_systems =
                    generator->get_mapped_transition_systems();
            merge_linear_transition_systems.insert(merge_linear_transition_systems.end(),
                                             mapped_transition_systems.begin(),
                                             mapped_transition_systems.end());
        }

        // If merging for least/most number of overall affected abstactions,
        // also include the non-mapped, i.e. internally affected transition systems
        // (always as to be linearly merged transition systems)
        if (external_merging == MERGE_FOR_ATOMIC) {
            const vector<int> &internally_affected_transition_systems =
                    generator->get_internally_affected_transition_systems();
            merge_linear_transition_systems.insert(merge_linear_transition_systems.end(),
                                             internally_affected_transition_systems.begin(),
                                             internally_affected_transition_systems.end());
        }

        // compute a merge tree
        assert(merge_order.empty());
        int number_of_transition_systems = transition_systems.size();
        int number_of_merges = 0;
        vector<int> merge_linear_indices;
        for (size_t cycle_no = 0; cycle_no < cycles.size(); ++cycle_no) {
            // go over the cycles and compute a non-linear merge order.
            const vector<int> &cycle = cycles[cycle_no];
            size_t abs_index_1 = cycle[0];
            for (size_t i = 1; i < cycle.size(); ++i) {
                size_t abs_index_2 = cycle[i];
                merge_order.push_back(make_pair(abs_index_1, abs_index_2));
                abs_index_1 = number_of_transition_systems + number_of_merges;
                ++number_of_merges;
            }
            if (external_merging == MERGE_FOR_ATOMIC) {
                // number_of_transition_systems + number_of_merges always is the *next*
                // position where a new merged transition system will be stored at.
                // here, we need the *last* position where the transition system
                // resulting from merging the cycle was stored, hence the -1.
                merge_linear_indices.push_back(number_of_transition_systems + number_of_merges - 1);
            }
        }

        if (external_merging == MERGE_FOR_ATOMIC) {
            // merge_linear_indices possibly contains transition systems that have been
            // non-linearly merged from information about cycles.
            // here we add transition systems that need to be merged linearly anyways
            merge_linear_indices.insert(merge_linear_indices.end(),
                                        merge_linear_transition_systems.begin(),
                                        merge_linear_transition_systems.end());

            // go over all transition systems that (now) need to be merged linearly
            size_t abs_index_1 = merge_linear_indices[0];
            for (size_t i = 1; i < merge_linear_indices.size(); ++i) {
                size_t abs_index_2 = merge_linear_indices[i];
                merge_order.push_back(make_pair(abs_index_1, abs_index_2));
                abs_index_1 = number_of_transition_systems + number_of_merges;
                ++number_of_merges;
            }
        }

        cout << "current number of transition systems " << number_of_transition_systems << endl;
        cout << "chosen internal merge order: " << endl;
        for (size_t i = 0; i < merge_order.size(); ++i) {
            cout << merge_order[i].first << ", " << merge_order[i].second << endl;
        }
    }

    return applied_symmetries;
}

void Symmetries::apply_symmetries(const vector<TransitionSystem *> &transition_systems,
                                  const vector<int> &generator_indices) const {
    if (get_num_generators() == 0) {
        cerr << "You first have to find symmetries before you can apply one of them!" << endl;
        exit_with(EXIT_CRITICAL_ERROR);
    }
    cout << "Creating equivalence relations from symmetries. [t=" << g_timer << "]" << endl;

    // Abstracting by the generators as follows:
    // Creating a graph with nodes being the abstract states and the edges represent the connections as given by the generators.
    // It seems like this process should better be done in all transition systems in parallel,
    // since we could exploit all the compositions of these transition systems this way as well.
    // Later we can reduce the first part - the transition system vars and save some space/time
    int num_states = get_sym_gen_info().num_abs_and_states;
    int num_transition_systems = get_sym_gen_info().num_transition_systems;
    // The graph is represented by vector of to_nodes for each node. (Change to sets?)
    vector<vector<int> > graph(num_states, vector<int>());
    for (int index = num_transition_systems; index < get_sym_gen_info().num_abs_and_states; ++index) {
        int abs_index = get_sym_gen_info().get_var_by_index(index);
        if (transition_systems[abs_index]) {
            for (size_t i = 0; i < generator_indices.size(); ++i) {
                // Going over the generators, for each just add the edges.
                if (get_symmetry_generator(generator_indices[i])->get_value(abs_index) == abs_index) {
                    // we only add an edge if the corresponding states belong to
                    // the same transition systems. in other words, we do not compute
                    // equivalence relations for mappings of transition systems, as
                    // these are not applied anyways.
                    int to_index = get_symmetry_generator(generator_indices[i])->get_value(index);
                    if (index != to_index)
                        graph[index].push_back(to_index);
                }
            }
        }
    }
    SCC scc(graph);
    const vector<vector<int> > &result = scc.get_result();

    // Generate final result. Going over the result, putting the nodes to their respective places.
    vector<vector<__gnu_cxx::slist<AbstractStateRef> >> equivalence_relations(transition_systems.size());
    for (size_t abs_index = 0; abs_index < transition_systems.size(); ++abs_index) {
        if (transition_systems[abs_index])
            equivalence_relations[abs_index].resize(result.size());
    }
    for (size_t eqiv=0; eqiv < result.size(); eqiv++) {
        for (size_t i=0; i < result[eqiv].size(); i++) {
            int idx = result[eqiv][i];
            if (idx < get_sym_gen_info().num_transition_systems)
                continue;
            pair<int, AbstractStateRef> vals = get_sym_gen_info().get_var_val_by_index(idx);
            equivalence_relations[vals.first][eqiv].push_front(vals.second);
        }
    }
    // Then, going over the outcome, removing empty equivalence classes.
    for (int eqiv=result.size(); eqiv > 0; --eqiv) {
        for (size_t i = 0; i < transition_systems.size(); ++i) {
            if (transition_systems[i] == 0)  //In case the transition system is empty
                continue;

            if (equivalence_relations[i][eqiv - 1].size() == 0)
                equivalence_relations[i].erase(equivalence_relations[i].begin() + eqiv - 1);
        }
    }

    cout << "==========================================================================================" << endl;
    cout << "Abstracting everything by the equivalence relations. [t=" << g_timer << "]" << endl;

    for (size_t i = 0; i < transition_systems.size(); ++i) {
        if (transition_systems[i] == 0)  //In case the transition system is empty
            continue;

        // Abstracting by the computed equivalence relation
        int equivalence_relation_size = equivalence_relations[i].size();
        if (equivalence_relation_size > transition_systems[i]->get_size()) {
            cerr << "Something is seriously wrong here!!" << endl;
            exit_with(EXIT_CRITICAL_ERROR);
        }
        if (equivalence_relation_size == transition_systems[i]->get_size()) {
            cout << transition_systems[i]->tag() << " not abstracted due to symmetries." << endl;
            continue;
        }
//      cout << "Abstracting from " << transition_systems[i]->size() << " to " << equivalence_relations[i].size() << " states!" << endl;
        cout << "Abstracting " << transition_systems[i]->tag() << endl;
        /*transition_systems[i]->dump();
        transition_systems[i]->dump_state();
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
        transition_systems[i]->apply_abstraction(equivalence_relations[i]);
    }
    cout << "Done abstracting. [t=" << g_timer << "]" << endl;
    cout << "==========================================================================================" << endl;
}

const SymmetryGenerator* Symmetries::get_symmetry_generator(int ind) const {
    assert(ind >= 0 && ind < get_num_generators());
    return gc.get_symmetry_generators()[ind];
}
