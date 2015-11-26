#ifndef MERGE_AND_SHRINK_SYMMETRIES_MS_GRAPH_CREATOR_H
#define MERGE_AND_SHRINK_SYMMETRIES_MS_GRAPH_CREATOR_H

#include <memory>

namespace bliss {
    class Digraph;
}
class FactoredTransitionSystem;
class Options;
class SymmetryGenerator;
class SymmetryGeneratorInfo;
class SymmetryGroup;

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

    void create_bliss_directed_graph(std::shared_ptr<FactoredTransitionSystem> fts,
                                     bliss::Digraph &bliss_graph,
                                     SymmetryGeneratorInfo *symmetry_generator_info);
public:
    explicit MSGraphCreator(const Options &options);
    ~MSGraphCreator();

    double compute_symmetries(std::shared_ptr<FactoredTransitionSystem> fts,
                              SymmetryGroup *symmetry_group,
                              SymmetryGeneratorInfo *symmetry_generator_info);
};

#endif
