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
        ATOMIC,
        LOCAL,
        NONE
    };
    SymmetriesForShrinking symmetries_for_shrinking;

    enum SymmetriesForMerging {
        SMALLEST,
        LARGEST
    };
    SymmetriesForMerging symmetries_for_merging;

    // search for local symmetries if true, for general ones if false
    bool build_stabilized_pdg;

    // the following serves for statistics output
    int atomic_symmetries; // symmetries affecting one abstraction
    int binary_symmetries; // symmetries affecting two abstractions
    int other_symmetries; // symmetries affecting more than two abstractions

    bool find_symmetries(const std::vector<Abstraction *>& abstractions);
    void apply_symmetries(const std::vector<Abstraction *> &abstractions,
                          const std::vector<int> &generator_indices) const;

    // TODO: replace by permutations wrapper object
    const SymmetryGenerator* get_symmetry_generator(int ind) const;
    const SymmetryGeneratorInfo &get_sym_gen_info() const { return gc.get_symmetry_generator_info(); }
    int get_num_generators() const { return gc.get_symmetry_generators().size(); }
public:
    explicit Symmetries(const Options &options);
    ~Symmetries() {}

    std::pair<int, int> find_and_apply_symmetries(std::vector<Abstraction *> &abstractions,
                                                  std::vector<int> &abs_to_merge);
    int get_atomic_symmetries() const {return atomic_symmetries; }
    int get_binary_symmetries() const {return binary_symmetries; }
    int get_other_symmetries() const {return other_symmetries; }
};

#endif
