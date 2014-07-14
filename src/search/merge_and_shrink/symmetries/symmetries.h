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
        LEAST_OVERALL_AFFECTED,
        MOST_OVERALL_AFFECTED,
        LEAST_MAPPED,
        MOST_MAPPED
    };
    SymmetriesForMerging symmetries_for_merging;

    enum InternalMerging {
        LINEAR,
        NON_LINEAR,
        NON_LINEAR_INCOMPLETE
    };
    InternalMerging internal_merging;

    // search for local symmetries if true, for general ones if false
    bool build_stabilized_pdg;

    void find_symmetries(const std::vector<Abstraction *>& abstractions);
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
};

#endif
