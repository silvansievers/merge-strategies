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
    // TODO: get rid of version or introduce an implementation for version 0
    int version;

    bool is_atomic_generator(const std::vector<Abstraction *> abstractions, int gen_index) const;
    bool find_atomic_symmetries(const std::vector<Abstraction *>& abstractions,
                                std::vector<std::vector<int> > &atomic_symmetries_by_affected_abs);

    bool find_symmetries(const std::vector<Abstraction *>& abstractions,
                         std::vector<std::set<int> > &affected_abstractions_by_generator,
                         std::vector<int> &atomic_generators);
    // TODO: there is a lot of duplicated code in here!
    void apply_symmetry(const std::vector<Abstraction *> &abstractions,
                        int generator_index) const;
    void apply_symmetries(const std::vector<Abstraction *> &abstractions,
                          const std::vector<int> &generator_indices) const;

    // TODO: replace by permutations wrapper object
    const Permutation* get_generator(int ind) const;
    const PermutationsWrapper &get_pw() const { return gc.get_permutations_wrapper(); }
    int get_num_generators() const { return gc.get_generators().size(); }
public:
    explicit Symmetries(const Options &options);
    ~Symmetries() {}

    bool find_and_apply_symmetries(const std::vector<Abstraction *> &abstractions,
                                        std::set<int> &abs_to_merge);
};

#endif
