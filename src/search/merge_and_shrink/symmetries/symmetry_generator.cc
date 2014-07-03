#include "symmetry_generator.h"

#include "../../globals.h"
#include "../../utilities.h"

#include <algorithm>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

using namespace std;

SymmetryGeneratorInfo::SymmetryGeneratorInfo() {
    reset();
}

void SymmetryGeneratorInfo::reset() {
    num_abstractions = -1;
    num_abs_and_states = -1;
    length = -1;
    var_by_val.clear();
    dom_sum_by_var.clear();
}

bool SymmetryGeneratorInfo::initialized() const {
    return num_abstractions != -1
            && num_abs_and_states != -1
            && length != -1
            && !var_by_val.empty()
            && !dom_sum_by_var.empty();
}

int SymmetryGeneratorInfo::get_var_by_index(const unsigned int ind) const {
    assert(initialized());
    if (ind < num_abstractions) {
        cout << "=====> WARNING!!!! Check that this is done on purpose!" << endl;
        return ind;
    }
    return var_by_val[ind-num_abstractions];
}

pair<int, AbstractStateRef> SymmetryGeneratorInfo::get_var_val_by_index(unsigned int ind) const {
    assert(initialized());
    if (ind < num_abstractions) {
        cout << "=====> Error!!!! index too low, in the variable part!" << endl;
        exit_with(EXIT_CRITICAL_ERROR);
    }

    int var =  var_by_val[ind-num_abstractions];
    int val = ind - dom_sum_by_var[var];

    if (val < 0) {
        cout << "=====> Error!!!! Problem with the index" << endl;
        exit_with(EXIT_CRITICAL_ERROR);
    }
//  cout << "=====================>" << var << " = " << val << endl;

    return make_pair(var, val);
}

unsigned int SymmetryGeneratorInfo::get_index_by_var_val_pair(int var, AbstractStateRef val) const {
    assert(initialized());
    return dom_sum_by_var[var] + val;
}

void SymmetryGeneratorInfo::dump() const {
    cout << "num abstractions: " << num_abstractions << endl;
    cout << "num abs and states: " << num_abs_and_states << endl;
    cout << "length: " << length << endl;
    cout << "var by val" << endl;
    cout << var_by_val << endl;
    cout << "dom sum by var" << endl;
    cout << dom_sum_by_var << endl;
}

void SymmetryGeneratorInfo::dump_var_by_val() const {
    int size = num_abs_and_states - num_abstractions;
    for (size_t i = 0; i < size; ++i) {
        cout << i << ": " << var_by_val[i];
        if (i != size - 1)
            cout << ", ";
    }
    cout << endl;
}




SymmetryGenerator::SymmetryGenerator(const SymmetryGeneratorInfo &sym_gen_info_,
                                     const unsigned int* automorphism,
                                     bool abstraction_stabilized_symmetry)
    : sym_gen_info(sym_gen_info_),
      identity_generator(true)/*,
      largest_cycle_size(0),
      largest_cycle_index(0)*/ {
    _allocate();

    int num_abstractions = sym_gen_info.num_abstractions;
    internally_affected.resize(num_abstractions, false);
    mapped.resize(num_abstractions, false);
    overall_affected.resize(num_abstractions, false);
    for (unsigned int from_index = 0; from_index < sym_gen_info.length; from_index++){
        if (from_index > sym_gen_info.num_abs_and_states) {
            cerr << "Symmetry generator index out of range" << endl;
            exit_with(EXIT_CRITICAL_ERROR);
        }

        int to_index = automorphism[from_index];
        value[from_index] = to_index;

        if (from_index != to_index) {
            identity_generator = false;
            if (from_index < num_abstractions) {
                // abstraction is mapped
                assert(to_index < num_abstractions);
                if (!mapped[from_index]) {
                    mapped[from_index] = true;
                    mapped_abstractions.push_back(from_index);
                }
                if (mapped[from_index] && internally_affected[from_index]) {
                    cerr << "Abstraction " << from_index << "both internally "
                         << "affected and mapped to another abstraction" << endl;
                    exit_with(EXIT_CRITICAL_ERROR);
                }
                if (!overall_affected[from_index]) {
                    overall_affected[from_index] = true;
                    overall_affected_abstractions.push_back(from_index);
                }
            } else {
                int from_abs_index = sym_gen_info.get_var_by_index(from_index);
                int to_abs_index = sym_gen_info_.get_var_by_index(to_index);
                if (!overall_affected[from_abs_index]) {
                    overall_affected[from_abs_index] = true;
                    overall_affected_abstractions.push_back(from_abs_index);
                }
                if (from_abs_index == to_abs_index) {
                    // abstraction affected internally
                    if (!internally_affected[from_abs_index]) {
                        internally_affected_abstractions.push_back(from_abs_index);
                        internally_affected[from_abs_index] = true;
                    }
                    if (mapped[from_abs_index] && internally_affected[from_abs_index]) {
                        cerr << "Abstraction " << from_abs_index << "both internally "
                             << "affected and mapped to another abstraction" << endl;
                        exit_with(EXIT_CRITICAL_ERROR);
                    }
                } else {
                    if (automorphism[from_abs_index] != to_abs_index) {
                        cerr << "State of abstraction mapped to state of another"
                             << " abstraction which differs from the abstractions'"
                             << " nodes mapping." << endl;
                        exit_with(EXIT_CRITICAL_ERROR);
                    }
                }
            }
        }
    }

    sort(internally_affected_abstractions.begin(), internally_affected_abstractions.end());
    sort(mapped_abstractions.begin(), mapped_abstractions.end());

    if (!abstraction_stabilized_symmetry)
        compute_cycles();
}

SymmetryGenerator::~SymmetryGenerator(){
    _deallocate();
}

void SymmetryGenerator::_allocate() {
    borrowed_buffer = false;
    value = new unsigned int[sym_gen_info.length];

    //reset_affected();
}

void SymmetryGenerator::_deallocate() {
    if (!borrowed_buffer) {
        delete[] value;
    }
}

void SymmetryGenerator::compute_cycles() {
    return;
    int num_abstractions = sym_gen_info.num_abstractions;
    vector<bool> marked(num_abstractions, false);
    for (size_t abs_index = 0; abs_index < num_abstractions; ++abs_index) {
        if (mapped[abs_index] && !marked[abs_index]) {
            marked[abs_index] = true;
            unsigned int to_index = get_value(abs_index);
            assert(to_index != abs_index);
            int from_index = abs_index;
            vector<int> cycle;
            cycle.push_back(from_index);
            while (to_index != abs_index) {
                marked[to_index] = true;
                cycle.push_back(to_index);
                from_index = to_index;
                to_index = get_value(from_index);
            }
            cycles.push_back(cycle);
            /*if (cycle.size > largest_cycle_size) {
                largest_cycle_size = cycle.size();
                largest_cycle_index = cycles.size() - 1;
            }*/
        }
    }
}

void SymmetryGenerator::get_mappings_for_cycles(vector<vector<pair<int, vector<int> > > > &mapping) const {
    mapping.reserve(cycles.size());
    for (size_t cycle_no = 0; cycle_no < cycles.size(); ++cycle_no) {
        const vector<int> &cycle = cycles[cycle_no];
        vector<pair<int, vector<int> > > cycle_mappings;
        cycle_mappings.reserve(cycle.size());
        cout << "cycle " << cycle_no << endl;
        cout << sym_gen_info.dom_sum_by_var << endl;
        dump_value();
        sym_gen_info.dump_var_by_val();
        for (size_t i = 0; i < cycle.size(); ++i) {
            size_t from_abs_index = cycle[i];
            size_t to_abs_index;
            if (i != cycle.size() - 1)
                to_abs_index = cycle[i + 1];
            else
                to_abs_index = cycle[0];
            cout << "abstraction " << from_abs_index << " -> " << to_abs_index << endl;
            vector<int> internal_abs_mapping;
            size_t value_index = 0;
            while (sym_gen_info.var_by_val[value_index] != from_abs_index) {
                // find starting index in value[] for abstraction from_abs_index
                ++value_index;
            }
            while (value_index < sym_gen_info.num_abs_and_states - sym_gen_info.num_abstractions && sym_gen_info.var_by_val[value_index] == from_abs_index) {
                // the entry x in var_by_val corresponds to the index
                // x + num_abstractions in value[]
                //cout << "value_index " << value_index << endl;
                int from_index = value_index + sym_gen_info.num_abstractions;
                //cout << "from_index " << from_index << endl;
                int to_index = value[from_index];
                //cout << "to_index " << to_index << endl;
                //assert(sym_gen_info.get_var_by_index(to_index) == to_abs_index);
                int to_abs_index = sym_gen_info.get_var_by_index(to_index);
                //cout << "to_abs_index " << to_abs_index << endl;
                if (i != cycle.size() - 1)
                    assert(to_abs_index == cycle[i + 1]);
                else
                    assert(to_abs_index == cycle[0]);
                int to_state_index = to_index - sym_gen_info.dom_sum_by_var[to_abs_index];
                cout << from_index - sym_gen_info.dom_sum_by_var[from_abs_index] << " -> " << to_state_index << endl;
                internal_abs_mapping.push_back(to_state_index);
                ++value_index;
            }
            cycle_mappings.push_back(make_pair(from_abs_index, internal_abs_mapping));
        }
        mapping.push_back(cycle_mappings);
    }
}

/*void Permutation::set_affected(unsigned int ind, unsigned int val) {
    if (ind < pw.num_abstractions || ind == val || ind >= pw.num_abs_and_states)
        return;

    int var = pw.get_var_by_index(ind);
    int to_var = pw.get_var_by_index(val);

//  cout << "Setting the affected variables - variable " << var << " is mapped to the variable " << to_var << endl;

    if (!affected[var]) {
        vars_affected.push_back(var);
        affected[var] = true;
    }
    if (!affected[to_var]) {
        vars_affected.push_back(to_var);
        affected[to_var] = true;
    }
    // Keeping the orig. var for each var.
    from_vars[to_var] = var;
}

void Permutation::reset_affected() {
    affected.assign(pw.num_abstractions, false);
    vars_affected.clear();
    from_vars.assign(pw.num_abstractions, -1);

    //affected_vars_cycles.clear();
}

void Permutation::finalize(){
    // Sorting the vector of affected variables
    ::sort(vars_affected.begin(), vars_affected.end());

    // Going over the vector from_vars of the mappings of the variables and finding cycles
//  affected_vars_cycles.clear();
    vector<bool> marked;
    marked.assign(pw.num_abstractions, false);
//  for (int i = 0; i < from_vars.size(); i++) {
//      cout << from_vars[i] << " ";
//  }
//  cout << endl;
    for (int i = 0; i < from_vars.size(); i++) {
        if (marked[i] || from_vars[i] == -1)
            continue;
//      if (from_vars[i] == i)
//          continue;

        int current = i;
        marked[current] = true;
        vector<int> cycle;
        cycle.push_back(current);

        while (from_vars[current] != i){
            current = from_vars[current];
            marked[current] = true;
            cycle.insert(cycle.begin(),current);
        }
        // Get here when from_vars[current] == i.
        affected_vars_cycles.push_back(cycle);
    }

    set_maximal_variables_cycle_size();
}*/

/*void Permutation::set_maximal_variables_cycle_size() {
    max_var_cycle_size = 0;
    for (int i=0; i < affected_vars_cycles.size(); i++) {
        if (max_var_cycle_size < affected_vars_cycles[i].size())
            max_var_cycle_size = affected_vars_cycles[i].size();
    }
}*/

bool SymmetryGenerator::identity() const{
    if (identity_generator)
        assert(internally_affected_abstractions.empty());
    return identity_generator;
}

unsigned int SymmetryGenerator::get_value(unsigned int ind) const {
//  check_index_range(ind);
    return value[ind];
}

/*
unsigned int Permutation::get_inverse_value(unsigned int ind) const {
//  check_index_range(ind);
    return inverse_value[ind];
}
*/

/*void Permutation::print_variables_by_cycles() const {

    cout << "Affected variables by cycles: " << endl;
    for (int i=0; i < affected_vars_cycles.size(); i++) {
        cout << "( " ;
        for (int j=0; j < affected_vars_cycles[i].size(); j++) {
            cout << affected_vars_cycles[i][j] << " ";
        }
        cout << ")  ";
    }
    cout << endl;
}*/

/*int Permutation::get_maximal_variables_cycle_size() const {
    return max_var_cycle_size;
}*/

/*int Permutation::calculate_number_variables_to_merge(bool linear_merge) const {
    int num_vars = 0;
    for (int i=0; i < affected_vars_cycles.size(); i++) {
        if (affected_vars_cycles[i].size() < 2)
            continue;

        num_vars += affected_vars_cycles[i].size();
        if (!linear_merge)
            num_vars--;
    }
    return num_vars;
}*/

/*void Permutation::print_cycle_notation() const {
//TODO revise
    vector<int> done;
    for (unsigned int i = pw.num_abstractions; i < pw.num_abs_and_states; i++){
        if (find(done.begin(), done.end(), i) == done.end()){
            unsigned int current = i;
            if(get_value(i) == i) continue; //don't print cycles of size 1

            pair<int, AbstractStateRef> varval = pw.get_var_val_by_index(i);
//          cout<<"("<< varval.first << "=" << (int) varval.second <<" ";
            cout<<"("<< g_fact_names[varval.first][(int) varval.second]  <<" ";

//          cout<<"("<<i<<" ";

            while(get_value(current) != i){
                done.push_back(current);
                current = get_value(current);

                pair<int, AbstractStateRef> currvarval = pw.get_var_val_by_index(current);
                // Silvan:
                //cout << currvarval.first << " " << currvarval.second << endl;
                cout<< g_fact_names[currvarval.first][(int) currvarval.second] <<" ";
//              cout<< currvarval.first << "=" << (int) currvarval.second <<" ";
//              cout<<current<<" ";

            }
            done.push_back(current);
            cout<<") ";
        }
    }


//    cout << endl << "Variables:  ";
//    for(int i = 0; i < vars_affected.size(); i++) cout << vars_affected[i] << "  ";
//    cout << endl << "Variables permuted:  ";

//    for(int i = 0; i < vars_affected.size(); i++) cout << from_vars[vars_affected[i]] << " -> " << vars_affected[i] << "  ";

//    cout << endl;
}*/

/*
string Permutation::get_cycle_notation() const {
    std::stringstream ss;
    ss << "  ";
    vector<int> done;
    for (unsigned int i = pw.num_abstractions; i < pw.num_abs_and_states; i++){
        if (find(done.begin(), done.end(), i) == done.end()){
            unsigned int current = i;
            if(get_value(i) == i) continue; //don't print cycles of size 1

            pair<int, AbstractStateRef> varval = get_var_val_by_index(i);
            ss << "( " << g_fact_names[varval.first][varval.second];
//          ss << "( " << i;
            while(get_value(current) != i){
                done.push_back(current);
                current = get_value(current);
                varval = get_var_val_by_index(current);
                ss << "  " << g_fact_names[varval.first][varval.second];
//              ss << "  " << current;
            }
            done.push_back(current);
            ss << " ), ";
        }
    }
    string first = ss.str();
    return first.substr(0, first.size()-2);
}
*/

void SymmetryGenerator::dump() const {
    for(unsigned int i = 0; i < sym_gen_info.length; i++){
        if (get_value(i) != i)
            cout << setw(4) << i;
    }
    cout << endl;
    for(unsigned int i = 0; i < sym_gen_info.length; i++){
        if (get_value(i) != i)
            cout << setw(4) << get_value(i);
    }
    cout << endl;
}

void SymmetryGenerator::dump_value() const {
    for (size_t i = 0; i < sym_gen_info.length; ++i) {
        cout << i << " -> " << value[i];
        if (i != sym_gen_info.length - 1)
            cout << ", ";
    }
    cout << endl;
}

void SymmetryGenerator::dump_all() const {
    cout << "values:" << endl;
    for(unsigned int i = 0; i < sym_gen_info.length; i++){
        cout << value[i] << ", ";
    }
    cout << endl;
    //cout << "vars affected" << endl;
    //cout << vars_affected << endl;
    //cout << "affected" << endl;
    //cout << affected << endl;
    cout << "borrowed buffer: " << borrowed_buffer << endl;
    cout << "identiy perm: " << identity_generator << endl;
    //cout << "from vars" << endl;
    //cout << from_vars << endl;
    /*cout << "affected vars cycles" << endl;
    for (size_t i = 0; i < affected_vars_cycles.size(); ++i) {
        cout << i << endl;
        cout << affected_vars_cycles[i] << endl;
    }
    cout << "max var cycle size: " << max_var_cycle_size << endl;*/
    sym_gen_info.dump();
}
