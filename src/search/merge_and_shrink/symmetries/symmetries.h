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
                         std::vector<std::set<int> > &non_trivially_affected_abstractions,
                         std::vector<bool> &atomic_symmetries,
                         std::vector<std::vector<int> > &atomic_symmetries_by_affected_abs,
                         bool find_atomic_symmetry);
    // TODO: there is a lot of duplicated code in here!
    void apply_symmetry(const std::vector<Abstraction *> &abstractions, int generator_index) const;
    void apply_symmetries(const std::vector<Abstraction *> &abstractions, const std::vector<int> &indices) const;

    void print_generators_stat() const;
    // TODO: replace by permutations wrapper object
    const Permutation* get_generator(int ind) const;
    const PermutationsWrapper &get_permutations_wrapper() const { return gc.get_permutations_wrapper(); }
public:
    explicit Symmetries(const Options &options);
    ~Symmetries() {}

    void find_and_apply_atomic_symmetries(const std::vector<Abstraction *> &abstractions);
    bool find_to_be_merged_abstractions(const std::vector<Abstraction *> &abstractions,
                                        std::set<int> &abs_to_merge);

    // TODO: make private
    int get_num_generators() const { return gc.get_generators().size(); }
};

#endif
