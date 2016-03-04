#ifndef MERGE_AND_SHRINK_SYMMETRIES_SYMMETRY_GROUP_H
#define MERGE_AND_SHRINK_SYMMETRIES_SYMMETRY_GROUP_H

#include <memory>
#include <vector>

namespace options {
class Options;
}

namespace merge_and_shrink {
class FactoredTransitionSystem;
class MSGraphCreator;
class SymmetryGenerator;
class SymmetryGeneratorInfo;

class SymmetryGroup {
    MSGraphCreator *gc;
    SymmetryGeneratorInfo *symmetry_generator_info;
    std::vector<const SymmetryGenerator*> symmetry_generators;

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

    void apply_symmetries(FactoredTransitionSystem &fts,
                          const std::vector<int> &generator_indices) const;
public:
    explicit SymmetryGroup(const options::Options &options);
    ~SymmetryGroup();

    // method used by add_automorphism
    void create_symmetry_generator(const unsigned int *automorphism);
    bool find_and_apply_symmetries(FactoredTransitionSystem &fts,
                                   std::vector<std::pair<int, int> > &merge_order);
    bool is_bliss_limit_reached() const {
        return bliss_limit_reached;
    }
    int get_num_generators() const {
        return symmetry_generators.size();
    }
    double get_bliss_time() const {
        return bliss_time;
    }
    void set_bliss_limit_reached() {
        bliss_limit_reached = true;
    }
};
}

#endif
