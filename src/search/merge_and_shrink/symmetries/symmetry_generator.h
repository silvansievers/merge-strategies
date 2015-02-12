#ifndef MERGE_AND_SHRINK_SYMMETRIES_SYMMETRY_GENERATOR_H
#define MERGE_AND_SHRINK_SYMMETRIES_SYMMETRY_GENERATOR_H

#include <vector>

struct SymmetryGeneratorInfo {
    int num_transition_systems;
    int num_abs_and_states;
    // Silvan: this vector ranges over the number of abstract states of this generator.
    // for every node it contains the number of the transition system this node
    // belongs to.
    std::vector<int> var_by_val;
    // Silvan:
    // index 0 contains the number of transition systems, i.e. the starting point in this
    // vector where the states of transition system 1 start.
    // index 1 contains the highest number (i.e. the last node) that is related to transition system 1.
    // index 2 etc.
    // size: number of transition systems
    std::vector<int> dom_sum_by_var;

    SymmetryGeneratorInfo();
    void reset();
    bool initialized() const;
    // Returns the transition system variable corresponding to the abstract state
    int get_var_by_index(const int val) const;
    std::pair<int, int> get_var_val_by_index(const int ind) const;
    int get_index_by_var_val_pair(const int var, const int val) const;
    void dump() const;
    void dump_var_by_val() const;
};

class SymmetryGenerator {
    const SymmetryGeneratorInfo *sym_gen_info;
    int* value;

    bool borrowed_buffer;
    bool identity_generator;

    void _allocate();
    void _deallocate();

    std::vector<bool> internally_affected;
    std::vector<int> internally_affected_transition_systems;
    std::vector<bool> mapped;
    std::vector<int> mapped_transition_systems;
    std::vector<bool> overall_affected;
    std::vector<int> overall_affected_transition_systems;
    //std::vector<std::vector<int> > cycles;

    //void compute_cycles();
public:
    SymmetryGenerator(const SymmetryGeneratorInfo *sym_gen_info,
                      const unsigned int* automorphism);
    ~SymmetryGenerator();

    bool identity() const;
    int get_value(int ind) const;

    const std::vector<int> &get_internally_affected_transition_systems() const {
        return internally_affected_transition_systems;
    }
    const std::vector<int> &get_mapped_transition_systems() const {
        return mapped_transition_systems;
    }
    const std::vector<int> &get_overall_affected_transition_systems() const {
        return overall_affected_transition_systems;
    }
    //const std::vector<std::vector<int> > &get_cycles() const {
    //    return cycles;
    //}
    void compute_cycles(std::vector<std::vector<int> > &cycles) const;
    //void get_mappings_for_cycles(std::vector<std::vector<std::pair<int, std::vector<int> > > > &mapping) const;
    bool internally_affects(int transition_system_index) const {
        return internally_affected[transition_system_index];
    }
    bool maps(int transition_system_index) const {
        return mapped[transition_system_index];
    }
    void dump() const;
    void dump_value() const;
    void dump_all() const;
};

#endif
