#ifndef MERGE_AND_SHRINK_SYMMETRIES_SYMMETRY_GROUP_H
#define MERGE_AND_SHRINK_SYMMETRIES_SYMMETRY_GROUP_H

#include "graph_creator.h"

#include <vector>

class Labels;
class Options;
class SymmetryGenerator;
class SymmetryGeneratorInfo;
class TransitionSystem;

//void add_automorphism(void*, unsigned int, const unsigned int *automorphism);

class SymmetryGroup {
    GraphCreator *gc;
    SymmetryGeneratorInfo *symmetry_generator_info;
    std::vector<const SymmetryGenerator*> symmetry_generators;

    int num_identity_generators;
    //int stop_after_false_generated;
    bool bliss_limit_reached;
    bool stop_after_no_symmetries;

    enum SymmetriesForShrinking {
        NO_SHRINKING,
        ATOMIC,
        LOCAL
    };
    SymmetriesForShrinking symmetries_for_shrinking;

    enum SymmetriesForMerging {
        NO_MERGING,
        SMALLEST,
        LARGEST
    };
    SymmetriesForMerging symmetries_for_merging;

    enum ExternalMerging {
        MERGE_FOR_ATOMIC,
        MERGE_FOR_LOCAL
    };
    ExternalMerging external_merging;

    enum InternalMerging {
        LINEAR,
        NON_LINEAR
    };
    InternalMerging internal_merging;

    double bliss_time; // elapsed bliss time

    void apply_symmetries(const std::vector<TransitionSystem *> &transition_systems,
                          const std::vector<int> &generator_indices) const;
public:
    explicit SymmetryGroup(const Options &options);
    ~SymmetryGroup();

    // method used by add_automorphism
    void create_symmetry_generator(const unsigned int *automorphism);
    bool find_and_apply_symmetries(const std::vector<TransitionSystem *> &transition_systems,
                                   std::vector<std::pair<int, int> > &merge_order);
    bool is_bliss_limit_reached() const {
        return bliss_limit_reached;
    }
    int get_num_generators() const {
        return symmetry_generators.size();
    }
    int get_num_identity_generators() const {
        return num_identity_generators;
    }
    double get_bliss_time() const {
        return bliss_time;
    }
    void set_bliss_limit_reached() {
        bliss_limit_reached = true;
    }
};

#endif
