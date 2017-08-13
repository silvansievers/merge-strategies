#include "symmetry_group.h"

#include "ms_graph_creator.h"
#include "symmetry_generator.h"

#include "../factored_transition_system.h"
#include "../transition_system.h"
#include "../types.h"

#include "../../option_parser.h"

#include "../../algorithms/sccs.h"

#include "../../utils/system.h"
#include "../../utils/timer.h"

#include <cassert>
#include <iostream>
#include <limits>

using namespace std;

namespace merge_and_shrink {
SymmetryGroup::SymmetryGroup(const Options &options)
    : bliss_limit_reached(false),
      stop_after_no_symmetries(options.get<bool>("stop_after_no_symmetries")),
      symmetries_for_merging(SymmetriesForMerging(options.get_enum("symmetries_for_merging"))),
      external_merging(ExternalMerging(options.get_enum("external_merging"))),
      internal_merging(InternalMerging(options.get_enum("internal_merging"))),
      bliss_time(0) {
    gc = new MSGraphCreator(options);
    symmetry_generator_info = new SymmetryGeneratorInfo();
}

SymmetryGroup::~SymmetryGroup() {
    for (size_t i = 0; i < symmetry_generators.size(); i++){
        delete symmetry_generators[i];
    }
    delete symmetry_generator_info;
    delete gc;
}

void SymmetryGroup::create_symmetry_generator(const unsigned int *automorphism) {
    SymmetryGenerator* symmetry_generator = new SymmetryGenerator(symmetry_generator_info, automorphism);
    if (!symmetry_generator->identity()) {
        symmetry_generators.push_back(symmetry_generator);
    } else {
        delete symmetry_generator;
    }
}

void SymmetryGroup::find_symmetries(const FactoredTransitionSystem &fts,
                                              vector<pair<int, int> > &merge_order) {
    bliss_time = gc->compute_symmetries(fts, this, symmetry_generator_info);
    delete gc;
    gc = 0;
    if (stop_after_no_symmetries && symmetry_generators.empty()) {
        bliss_limit_reached = true;
    }
    if (symmetry_generators.empty() || bliss_limit_reached) {
        return;
    }


    int chosen_generator_for_merging = -1;
    int smallest_generator_affected_transition_systems_size = numeric_limits<int>::max();
    int smallest_generator_mapped_transition_systems_size = numeric_limits<int>::max();
    int largest_generator_affected_transition_systems_size = 0;
    int largest_generator_mapped_transition_systems_size = 0;

    /*
      Go over the generators and classify them into atomic, local or general
      ones. Also compute the generator that affects or maps the smallest or
      largest number of transition systems (depending on the chosen options).
    */
    for (int generator_index = 0; generator_index < get_num_generators(); ++generator_index) {
        const SymmetryGenerator *generator = symmetry_generators[generator_index];
        const vector<int> &mapped_transition_systems =
            generator->get_mapped_transition_systems();
        const vector<int> &overall_affected_transition_systems =
            generator->get_overall_affected_transition_systems();

        int number_overall_affected_transition_systems = overall_affected_transition_systems.size();
        if (number_overall_affected_transition_systems < 1) {
            cerr << "Something is wrong! The generator is the identity generator." << endl;
            utils::exit_with(utils::ExitCode::CRITICAL_ERROR);
        }
        if (number_overall_affected_transition_systems > 1) {
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
        if (number_mapped_transition_systems > 0) {
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

        // Dump generator properties -- warning! lots of output!
//        const vector<int> &internally_affected_transition_systems =
//            generator->get_internally_affected_transition_systems();
//        cout << "Generator " << generator_index << endl;
//        for (size_t i = 0; i < mapped_transition_systems.size(); ++i) {
//            int abs_index = mapped_transition_systems[i];
//            int to_index = generator->get_value(abs_index);
//            cout << transition_systems[abs_index]->tag() << " mapped to " <<
//                    transition_systems[to_index]->tag();
//            if (generator->internally_affects(abs_index))
//                cout << " (and also internally affected)";
//            cout << endl;
//        }
//        for (size_t i = 0; i < internally_affected_transition_systems.size(); ++i) {
//            int abs_index = internally_affected_transition_systems[i];
//            assert(!generator->maps(abs_index));
//            cout << transition_systems[abs_index]->tag() << " internally affected" << endl;
//        }
    }

    /*
      Determine a merge strategy as follows:
      If we merge towards local symmetries, we only need to merge all mapped
      transition systems. They are always merged according to the cycle(s)
      of those mapped transition systems. If we merge towards atomic symmetries,
      it depends on whether we want to merge all the affected transition systems
      linearly or non-linearly: if linearly, we put them in the order given by
      fast downward. If non-linearly, we first merge all cycles (as when merging
      for local symmetries), and then linearly the products of these merges
      and the remaining other (internally) affected transition systems, again
      as given by fast downward.
    */
    if (symmetries_for_merging != NO_MERGING && chosen_generator_for_merging != -1) {
        vector<vector<int> > cycles;
        vector<int> merge_linear_transition_systems;
        const SymmetryGenerator *generator = symmetry_generators[chosen_generator_for_merging];

        // Always include all mapped transition systems
        if (internal_merging == NON_LINEAR
                || external_merging == MERGE_FOR_LOCAL) {
            /*
              If the internal merge strategy is non linear or we only want
              to merge every cycle (non linearly), we need to
              compute the actual cycles of transition system mappings.
            */
            generator->compute_cycles(cycles);
        } else if (internal_merging == LINEAR) {
            /*
              If the internal merge strategy is linear, we simply collect
              all mapped transition systems (i.e. the same transition systems
              as above, but we do not compute cycle information)
            */
            const vector<int> &mapped_transition_systems =
                generator->get_mapped_transition_systems();
            merge_linear_transition_systems.insert(merge_linear_transition_systems.end(),
                                                   mapped_transition_systems.begin(),
                                                   mapped_transition_systems.end());
        }

        /*
          If merging for atomic symmetries, we need to include the internally
          affected transition systems (to be merge linearly after the mapped
          ones).
        */
        if (external_merging == MERGE_FOR_ATOMIC) {
            const vector<int> &internally_affected_transition_systems =
                 generator->get_internally_affected_transition_systems();
            merge_linear_transition_systems.insert(merge_linear_transition_systems.end(),
                                                   internally_affected_transition_systems.begin(),
                                                   internally_affected_transition_systems.end());
        }

        /*
          Here we compute the actual merge tree: if cycles is non-empty, we
          start by merging all mapped transition systems cycle-wise. At the
          same time, we need to remember the indices of the transition systems
          resulting from those merges in case we want to merge those afterwards
          (when merging for atomic symmetries).
        */
        assert(merge_order.empty());
        int number_of_transition_systems = fts.get_size();
        int number_of_merges = 0;
        vector<int> merge_linear_indices;
        for (size_t cycle_no = 0; cycle_no < cycles.size(); ++cycle_no) {
            const vector<int> &cycle = cycles[cycle_no];
            size_t abs_index_1 = cycle[0];
            for (size_t i = 1; i < cycle.size(); ++i) {
                size_t abs_index_2 = cycle[i];
                merge_order.push_back(make_pair(abs_index_1, abs_index_2));
                abs_index_1 = number_of_transition_systems + number_of_merges;
                ++number_of_merges;
            }
            if (external_merging == MERGE_FOR_ATOMIC) {
                /*
                  number_of_transition_systems + number_of_merges always is the *next*
                  position where a new merged transition system will be stored at.
                  here, we need the *last* position where the transition system
                  resulting from merging the cycle was stored, hence the -1.
                */
                merge_linear_indices.push_back(number_of_transition_systems + number_of_merges - 1);
            }
        }

        /*
          Here we include the transition systems that need to be merged linearly
          in the merge tree. Those are the internally affected ones if they
          need to be merged, and the products of the merged cycles if
          applicable (merge_linear_transition_systems).
        */
        if (external_merging == MERGE_FOR_ATOMIC) {
            merge_linear_indices.insert(merge_linear_indices.end(),
                                        merge_linear_transition_systems.begin(),
                                        merge_linear_transition_systems.end());

            size_t abs_index_1 = merge_linear_indices[0];
            for (size_t i = 1; i < merge_linear_indices.size(); ++i) {
                size_t abs_index_2 = merge_linear_indices[i];
                merge_order.push_back(make_pair(abs_index_1, abs_index_2));
                abs_index_1 = number_of_transition_systems + number_of_merges;
                ++number_of_merges;
            }
        }

        cout << "Chosen internal merge order: " << endl;
        for (size_t i = 0; i < merge_order.size(); ++i) {
            cout << merge_order[i].first << ", " << merge_order[i].second << endl;
        }
    }
}
}
