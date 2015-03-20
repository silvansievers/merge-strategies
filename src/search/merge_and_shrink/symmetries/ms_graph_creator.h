#ifndef MERGE_AND_SHRINK_SYMMETRIES_MS_GRAPH_CREATOR_H
#define MERGE_AND_SHRINK_SYMMETRIES_MS_GRAPH_CREATOR_H

#include <vector>

namespace bliss {
    class Digraph;
}
class Options;
class SymmetryGenerator;
class SymmetryGeneratorInfo;
class SymmetryGroup;
class TransitionSystem;

/**
 * This class is using bliss for finding symmetries of the given set of transition systems.
 */

class MSGraphCreator {
    enum color_t {
        TRANSITION_SYSTEM_VERTEX,
        ABSTRACT_STATE_VERTEX,
        GOAL_VERTEX,
        TRANSITION_VERTEX,
        LABEL_GROUP_VERTEX,
        LABEL_VERTEX,
        INITIAL_VERTEX
    };

    // Options
    bool debug; //generate dot-readable output
    bool stabilize_transition_systems;
    double bliss_time_limit;

    void create_bliss_directed_graph(const std::vector<TransitionSystem *>& transition_systems,
                                     bliss::Digraph &bliss_graph,
                                     SymmetryGeneratorInfo *symmetry_generator_info);
public:
    explicit MSGraphCreator(const Options &options);
    ~MSGraphCreator();

    double compute_symmetries(const std::vector<TransitionSystem *>& transition_systems,
                              SymmetryGroup *symmetry_group,
                              SymmetryGeneratorInfo *symmetry_generator_info);
};

#endif
