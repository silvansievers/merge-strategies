#ifndef MERGE_AND_SHRINK_SYMMETRIES_PERMUTATION_H
#define MERGE_AND_SHRINK_SYMMETRIES_PERMUTATION_H

typedef int AbstractStateRef; // TODO: duplicated from shrink_strategy.h

//#include <utility>
#include <vector>

struct PermutationsWrapper {
    int num_abstractions;
    unsigned int num_abs_and_states;
    unsigned int length;
    // Silvan: this vector ranges over the total number of indices of this permutation.
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

    PermutationsWrapper();
    void reset();
    bool initialized() const;
    // Returns the abstraction variable corresponding to the abstract state
    int get_var_by_index(const unsigned int val) const;
    std::pair<int, AbstractStateRef> get_var_val_by_index(const unsigned int ind) const;
    unsigned int get_index_by_var_val_pair(const int var, const AbstractStateRef val) const;
    void dump() const;
};

class Permutation {
    const PermutationsWrapper &pw;
    unsigned int* value;
    //unsigned int* inverse_value;
    //std::vector<int> vars_affected;
    //std::vector<bool> affected;
    bool borrowed_buffer;
    bool identity_perm;
    // Need to keep the connection between affected vars, ie which var goes into which.
    //std::vector<int> from_vars;
    // Affected vars by cycles
    //std::vector<std::vector<int> > affected_vars_cycles;
    //int max_var_cycle_size;

    void _allocate();
    void _deallocate();
    void set_value(unsigned int ind, unsigned int val);
    //void set_affected(unsigned int ind, unsigned int val);
    //void reset_affected();
    //void finalize();
    //void set_maximal_variables_cycle_size();
public:
    Permutation(const PermutationsWrapper &pw, const unsigned int* full_permutation);
    ~Permutation();

    bool identity() const;
    unsigned int get_value(unsigned int ind) const;
    //unsigned int get_inverse_value(unsigned int ind) const;
    //void print_variables_by_cycles() const;
    //int get_maximal_variables_cycle_size() const;
    //int calculate_number_variables_to_merge(bool linear_merge) const;
    //void print_cycle_notation() const;
    //string get_cycle_notation() const;
    void dump() const;
    void dump_all() const;
};

#endif // MERGE_AND_SHRINK_SYMMETRIES_PERMUTATION_H
