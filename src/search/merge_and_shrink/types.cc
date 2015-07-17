#include "types.h"

#include "../globals.h"

#include <utility>

using namespace std;
using namespace mst;

var_set_t mst::singleton(const var_t &var) {
    var_set_t singleton;
    singleton.insert(var);
    return singleton;
}

mutex_set_t mst::get_mutex_pairs() {
    mutex_set_t mutex_set;
    /* get mutex pairs */
    /* convert from the (first_variable, first_value) indexed representation */
    /* to the (first_variable, second_variable) indexed representation */
    for (size_t u = 0; u < g_inconsistent_facts.size(); ++u) {
        for (size_t value_u = 0; value_u < g_inconsistent_facts[u].size();
             ++value_u) {
            for (set<pair<int, int> >::iterator
                 i = g_inconsistent_facts[u][value_u].begin();
                 i != g_inconsistent_facts[u][value_u].end(); ++i) {
                size_t v = i->first;
                if (u == v)
                    continue;
                size_t value_v = i->second;

                var_set_t Variables;
                Variables.insert(u);
                Variables.insert(v);

                /* variables are sorted in increasing order (std::set defult)*/
                /* keep the values (std::vector) in the same order */
                /* the following process works only for mutex pairs */
                vector<value_t> Values;
                if (u < v) {
                    Values.push_back(value_u);
                    Values.push_back(value_v);
                } else {
                    Values.push_back(value_v);
                    Values.push_back(value_u);
                }
                /* if it is a new variable sets */
                if (!mutex_set.count(Variables)) {
                    set<vector<value_t> > vector_values;
                    vector_values.insert(Values);
                    mutex_set.insert(
                        pair<var_set_t, set<vector<value_t> > >(
                            Variables, vector_values));
                } else {
                    mutex_set.find(Variables)->second.insert(Values);
                }
            }
        }
    }
    return mutex_set;
}

set<var_set_t> mst::get_mutex_pairs_var() {
    set<var_set_t> mutex_var_set;
    /* get mutex pairs, keep the variables but not the values */
    for (size_t u = 0; u < g_inconsistent_facts.size(); ++u) {
        for (size_t value_u = 0; value_u < g_inconsistent_facts[u].size();
             ++value_u) {
            for (set<pair<int, int> >::iterator
                 i = g_inconsistent_facts[u][value_u].begin();
                 i != g_inconsistent_facts[u][value_u].end(); ++i) {
                size_t v = i->first;
                if (u == v)
                    continue;

                var_set_t Variables;
                Variables.insert(u);
                Variables.insert(v);

                mutex_var_set.insert(Variables);
            }
        }
    }
    return mutex_var_set;
}

var_relation_t mst::get_mutex_pairs_relation() {
    const size_t n = g_variable_domain.size();
    vector<var_t> row(n, 0);
    var_relation_t mutex_pair_relation(n, row);

    for (size_t u = 0; u < g_inconsistent_facts.size(); ++u) {
        for (size_t val_u = 0; val_u < g_inconsistent_facts[u].size();
             ++val_u) {
            for (set<pair<int, int> >::iterator
                 v_pair = g_inconsistent_facts[u][val_u].begin();
                 v_pair != g_inconsistent_facts[u][val_u].end();
                 ++v_pair) {
                var_set_t p;
                p.insert(u);
                p.insert(v_pair->first);
                mutex_pair_relation[u][v_pair->first]++;
                mutex_pair_relation[v_pair->first][u]++;
            }
        }
    }

    return mutex_pair_relation;
}
