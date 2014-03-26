#ifndef MERGE_AND_SHRINK_SYMMETRIES_GRAPH_CREATOR_H
#define MERGE_AND_SHRINK_SYMMETRIES_GRAPH_CREATOR_H

#include "permutation.h"

#include "../../bliss/graph.hh"

#include <vector>

class Abstraction;

void add_permutation(void*, unsigned int, const unsigned int *);

/**
 * This class is using bliss for finding symmetries of the given set of abstractions.
 */

class GraphCreator {
    enum color_t {
        ABSTRACTION_VERTEX,
        ABS_STATE_VERTEX,
        GOAL_VERTEX,
        LABEL_VERTEX,
        INITIAL_VERTEX
    };

    //int time_bound;
    //int generators_bound;
    bool debug; //generate dot-readable output
    int num_identity_generators;
    //int stop_after_false_generated;

    // TODO: some thouhths on this: maybe we can use the permutations_wrapper
    // to actually store all permutations. This would justify the name and
    // possibly simplify the handling of all permutations.
    std::vector<const Permutation*> generators; // the generators for the automorphism
    PermutationsWrapper permutations_wrapper;

    bliss::Digraph* create_bliss_graph(const std::vector<Abstraction *>& abstractions, bool stabilize_abstractions);

    void delete_generators();
public:
    GraphCreator(bool debug);
    ~GraphCreator();

    // method used by add_permutation
    void add_generator(const unsigned int *full_perm);

    void compute_generators(const std::vector<Abstraction *>& abstractions, bool stabilize_abstractions);
    const std::vector<const Permutation*>& get_generators () const { return generators; }
    const PermutationsWrapper &get_permutations_wrapper() const { return permutations_wrapper; }
};

#endif
