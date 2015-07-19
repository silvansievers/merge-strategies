#ifndef TYPES_H
#define TYPES_H

#include <set>
#include <map>
#include <vector>

#include <limits>

/** Merge-and-Shrink Types */
namespace mst {
typedef int var_t;
typedef int value_t;
typedef std::set<var_t> var_set_t;

var_set_t singleton(const var_t &var);

typedef std::map<var_set_t, std::set<std::vector<value_t> > > mutex_set_t;
mutex_set_t get_mutex_pairs();

std::set<var_set_t> get_mutex_pairs_var();

typedef std::vector<std::vector<int> > var_relation_t;
var_relation_t get_mutex_pairs_relation();
}

#endif // TYPES_H
