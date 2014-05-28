#ifndef MERGE_AND_SHRINK_SYMMETRIES_SYMMETRIES_H
#define MERGE_AND_SHRINK_SYMMETRIES_SYMMETRIES_H

#include "graph_creator.h"

#include "../abstraction.h"

#include <set>
#include <vector>

class Labels;
class Options;
class Permutation;
class PermutationsWrapper;

class Symmetries {
    // TODO: get rid of gc and have the permutations wrapper object instead?
    // This would store the actual generators as well.
    GraphCreator gc;

    enum TypeOfSymmetries {
        ATOMIC,
        LOCAL,
        ANY
    };
    TypeOfSymmetries type_of_symmetries;
    bool build_stabilized_pdg;

    // the following serves for statistics output
    int atomic_symmetries; // symmetries affecting one abstraction
    int binary_symmetries; // symmetries affecting two abstractions
    int other_symmetries; // symmetries affecting more than two abstractions

    bool find_symmetries(const std::vector<Abstraction *>& abstractions,
                         std::vector<std::set<int> > &affected_abstractions_by_generator,
                         std::vector<int> &atomic_generators,
                         std::vector<int> &local_generators);
    void apply_symmetries(const std::vector<Abstraction *> &abstractions,
                          const std::vector<int> &generator_indices) const;

    // TODO: replace by permutations wrapper object
    const Permutation* get_generator(int ind) const;
    const PermutationsWrapper &get_pw() const { return gc.get_permutations_wrapper(); }
    int get_num_generators() const { return gc.get_generators().size(); }
public:
    explicit Symmetries(const Options &options);
    ~Symmetries() {}

    int find_and_apply_symmetries(const std::vector<Abstraction *> &abstractions,
                                  std::set<int> &abs_to_merge);
    int get_atomic_symmetries() const {return atomic_symmetries; }
    int get_binary_symmetries() const {return binary_symmetries; }
    int get_other_symmetries() const {return other_symmetries; }
};

#endif
