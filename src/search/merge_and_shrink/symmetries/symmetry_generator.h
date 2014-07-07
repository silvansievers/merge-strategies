#ifndef MERGE_AND_SHRINK_SYMMETRIES_SYMMETRY_GENERATOR_H
#define MERGE_AND_SHRINK_SYMMETRIES_SYMMETRY_GENERATOR_H

typedef int AbstractStateRef; // TODO: duplicated from shrink_strategy.h

//#include <utility>
#include <set>
#include <vector>

struct SymmetryGeneratorInfo {
    int num_abstractions;
    unsigned int num_abs_and_states;
    unsigned int length;
    // Silvan: this vector ranges over the number of abstract states of this generator.
    // for every node it contains the number of the abstraction this node
    // belongs to.
    std::vector<int> var_by_val;
    // Silvan:
    // index 0 contains the number of abstractions, i.e. the starting point in this
    // vector where the states of abstraction 1 start.
    // index 1 contains the highest number (i.e. the last node) that is related to abstraction 1.
    // index 2 etc.
    // size: number of abstractions
    std::vector<unsigned int> dom_sum_by_var;

    SymmetryGeneratorInfo();
    void reset();
    bool initialized() const;
    // Returns the abstraction variable corresponding to the abstract state
    int get_var_by_index(const unsigned int val) const;
    std::pair<int, AbstractStateRef> get_var_val_by_index(const unsigned int ind) const;
    unsigned int get_index_by_var_val_pair(const int var, const AbstractStateRef val) const;
    void dump() const;
    void dump_var_by_val() const;
};

class SymmetryGenerator {
    const SymmetryGeneratorInfo &sym_gen_info;
    unsigned int* value;

    bool borrowed_buffer;
    bool identity_generator;

    void _allocate();
    void _deallocate();

    std::vector<bool> internally_affected;
    std::vector<int> internally_affected_abstractions;
    std::vector<bool> mapped;
    std::vector<int> mapped_abstractions;
    std::vector<bool> overall_affected;
    std::vector<int> overall_affected_abstractions;
    //std::vector<std::vector<int> > cycles;

    //void compute_cycles();
public:
    SymmetryGenerator(const SymmetryGeneratorInfo &sym_gen_info,
                      const unsigned int* automorphism,
                      bool abstraction_stabilized_symmetry);
    ~SymmetryGenerator();

    bool identity() const;
    unsigned int get_value(unsigned int ind) const;

    const std::vector<int> &get_internally_affected_abstractions() const {
        return internally_affected_abstractions;
    }
    const std::vector<int> &get_mapped_abstractions() const {
        return mapped_abstractions;
    }
    const std::vector<int> &get_overall_affected_abstractions() const {
        return overall_affected_abstractions;
    }
    //const std::vector<std::vector<int> > &get_cycles() const {
    //    return cycles;
    //}
    void compute_cycles(std::vector<std::vector<int> > cycles) const;
    //void get_mappings_for_cycles(std::vector<std::vector<std::pair<int, std::vector<int> > > > &mapping) const;
    bool internally_affects(int abs_index) const {return internally_affected[abs_index]; }
    bool maps(int abs_index) const {return mapped[abs_index]; }
    void dump() const;
    void dump_value() const;
    void dump_all() const;
};

#endif
