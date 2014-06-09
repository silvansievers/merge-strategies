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




SymmetryGenerator::SymmetryGenerator(const SymmetryGeneratorInfo &pw_, const unsigned int* symmetry_mapping)
    : pw(pw_), identity_generator(true)/*, max_var_cycle_size(-1)*/ {
//  cout << "Allocating" << endl;
    _allocate();
//  cout << "Setting values" << endl;
    for (unsigned int i = 0; i < pw.length; i++){
        set_value(i,symmetry_mapping[i]);
        if (i != symmetry_mapping[i])
            identity_generator = false;
    }
//  cout << "Finalizing" << endl;
    //finalize();
//  cout << "Done" << endl;
}

SymmetryGenerator::~SymmetryGenerator(){
    _deallocate();
}

void SymmetryGenerator::_allocate() {
    borrowed_buffer = false;
    value = new unsigned int[pw.length];

    //reset_affected();
}

void SymmetryGenerator::_deallocate() {
    if (!borrowed_buffer) {
        delete[] value;
    }
}

void SymmetryGenerator::set_value(unsigned int ind, unsigned int val) {
    value[ind] = val;
//  inverse_value[val] = ind;
    //set_affected(ind, val);
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
//  return vars_affected.size() == 0;
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
    for(unsigned int i = 0; i < pw.length; i++){
        if (get_value(i) != i)
            cout << setw(4) << i;
    }
    cout << endl;
    for(unsigned int i = 0; i < pw.length; i++){
        if (get_value(i) != i)
            cout << setw(4) << get_value(i);
    }
    cout << endl;
}

void SymmetryGenerator::dump_all() const {
    cout << "values:" << endl;
    for(unsigned int i = 0; i < pw.length; i++){
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
    pw.dump();
}
