#ifndef MERGE_AND_SHRINK_SYMMETRIES_GRAPH_CREATOR_H
#define MERGE_AND_SHRINK_SYMMETRIES_GRAPH_CREATOR_H

#include "symmetry_generator.h"

#include "../../bliss/graph.hh"

#include <vector>

class Options;
class TransitionSystem;

void add_automorphism(void*, unsigned int, const unsigned int *automorphism);

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

    // Options
    bool debug; //generate dot-readable output
    bool build_stabilized_pdg;
    double bliss_time_limit;
    bool stop_after_no_symmetries;

    int num_identity_generators;
    //int stop_after_false_generated;
    bool bliss_limit_reached;

    std::vector<const SymmetryGenerator*> symmetry_generators;
    SymmetryGeneratorInfo symmetry_generator_info;

    void create_bliss_graph(const std::vector<TransitionSystem *>& abstractions,
                            bliss::Digraph &bliss_graph);

    void delete_generators();
public:
    explicit GraphCreator(const Options &options);
    ~GraphCreator();

    // method used by add_automorphism
    void create_symmetry_generator(const unsigned int *automorphism);

    double compute_generators(const std::vector<TransitionSystem *>& abstractions);
    const std::vector<const SymmetryGenerator*>& get_symmetry_generators () const {
        return symmetry_generators;
    }
    const SymmetryGeneratorInfo &get_symmetry_generator_info() const {
        return symmetry_generator_info;
    }
    bool is_bliss_limit_reached() const {
        return bliss_limit_reached;
    }
};

#endif
