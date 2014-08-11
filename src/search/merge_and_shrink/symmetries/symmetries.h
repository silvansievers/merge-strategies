#ifndef MERGE_AND_SHRINK_SYMMETRIES_SYMMETRIES_H
#define MERGE_AND_SHRINK_SYMMETRIES_SYMMETRIES_H

#include "graph_creator.h"

#include "../abstraction.h"

#include <vector>

class Labels;
class Options;
class SymmetryGenerator;
class SymmetryGeneratorInfo;

class Symmetries {
    // TODO: get rid of gc and have the permutations wrapper object instead?
    // This would store the actual generators as well.
    GraphCreator gc;

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

    void apply_symmetries(const std::vector<Abstraction *> &abstractions,
                          const std::vector<int> &generator_indices) const;

    // TODO: replace by permutations wrapper object
    const SymmetryGenerator* get_symmetry_generator(int ind) const;
    const SymmetryGeneratorInfo &get_sym_gen_info() const { return gc.get_symmetry_generator_info(); }
    int get_num_generators() const { return gc.get_symmetry_generators().size(); }
public:
    explicit Symmetries(const Options &options);
    ~Symmetries() {}

    bool find_and_apply_symmetries(std::vector<Abstraction *> &abstractions,
                                   std::vector<std::pair<int, int> > &merge_order);
    bool is_bliss_limit_reached() const {return gc.is_bliss_limit_reached(); }
    double get_bliss_time() const {return bliss_time; }
};

#endif
