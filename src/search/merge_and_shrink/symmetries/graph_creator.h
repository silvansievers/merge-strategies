#ifndef MERGE_AND_SHRINK_SYMMETRIES_GRAPH_CREATOR_H
#define MERGE_AND_SHRINK_SYMMETRIES_GRAPH_CREATOR_H

#include "symmetry_generator.h"

#include "../../bliss/graph.hh"

#include <vector>

class Abstraction;
class Options;

void add_symmetry(void*, unsigned int, const unsigned int *);

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
    bool build_stabilized_pdg;
    int num_identity_generators;
    //int stop_after_false_generated;

    std::vector<const SymmetryGenerator*> symmetry_generators; // the generators for the automorphism
    SymmetryGeneratorInfo symmetry_generator_info;

    bliss::Digraph* create_bliss_graph(const std::vector<Abstraction *>& abstractions);

    void delete_generators();
public:
    explicit GraphCreator(const Options &options);
    ~GraphCreator();

    // method used by add_symmetry
    void add_symmetry_generator(const unsigned int *symmetry_mapping);

    void compute_generators(const std::vector<Abstraction *>& abstractions);
    const std::vector<const SymmetryGenerator*>& get_symmetry_generators () const { return symmetry_generators; }
    const SymmetryGeneratorInfo &get_symmetry_generator_info() const { return symmetry_generator_info; }
};

#endif
