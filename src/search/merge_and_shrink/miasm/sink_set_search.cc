#include "sink_set_search.h"

#include "miasm_mas.h"
#include "merge_miasm.h"

#include "../factored_transition_system.h"
#include "../transition_system.h"

#include "../../causal_graph.h"
#include "../../scc.h"

#include "get_rss.h"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include <algorithm>
#include <iostream>


using namespace std;

namespace MergeAndShrink {
using namespace mst;

SinkSetSearch::SinkSetSearch(const Options &opts, const shared_ptr<AbstractTask> task)
    : miasm_abstraction(
          opts.get<MiasmAbstraction *>(MiasmAbstraction::option_key())),
      task(task),
      task_proxy(*task),
      causal_graph(task_proxy.get_causal_graph()),
      time_limit(opts.get<double>(OptTimeLimit::opt_key())),
      memory_limit(opts.get<int>(OptMemoryLimit::opt_key())),
      size_limit(opts.get<int>(OptSizeLimit::opt_key())),
      clique_limit(opts.get<int>(OptCliqueLimit::opt_key())),
      opt_prior(opts.get_enum(EnumPriority::option_key())),
      opt_expa(opts.get_enum(EnumExpand::option_key())),
      opt_gain(opts.get_enum(EnumGain::option_key())),
      opt_prune(opts.get_enum(EnumPrune::option_key())),
      pq(ComparatorSTLPriorityQueue(task, &vsir, &opt_prior)) {
//    cerr << __PRETTY_FUNCTION__ << endl;
//    dump_options(cerr, "\n    ");
}

bool SinkSetSearch::time_limit_exceeded() {
    return timer() > time_limit;
}

bool SinkSetSearch::memory_limit_exceeded() {
    return getCurrentRSS() > memory_limit;
}

void SinkSetSearch::dump_options(ostream &os, const string sep) const {
    os << sep << OptTimeLimit::opt_key() << " = " << time_limit
       << sep << OptSizeLimit::opt_key() << " = " << size_limit
       << sep << OptCliqueLimit::opt_key() << " = " << clique_limit
       << sep << endl;
    os << sep << EnumPriority::option_key() << " = "
       << EnumPriority::C[opt_prior]
       << sep << EnumExpand::option_key() << " = "
       << EnumExpand::C[opt_expa]
       << sep << EnumGain::option_key() << " = "
       << EnumGain::C[opt_gain]
       << sep << EnumPrune::option_key() << " = "
       << EnumPrune::C[opt_prune]
       << sep << endl;
}

void SinkSetSearch::get_sink_set(vector<var_set_t> &sink_set) {
    std::sort(sink_set_idx.begin(), sink_set_idx.end(),
              ComparatorVarSet(task, &vsir, VarSetCmpType::BY_RATIO));

//    for (size_t i = 0; i < sink_set_idx.size(); i++) {
//        cerr << vsir[sink_set_idx[i]].variables << ": "
//             << vsir[sink_set_idx[i]].ratio << " "
//             << vsir[sink_set_idx[i]].gain << endl;
//    }

    for (size_t i = 0; i < sink_set_idx.size(); i++) {
        sink_set.push_back(vsir[sink_set_idx[i]].variables);
    }
}

void SinkSetSearch::reset() {
    vsir.clear();
    visited.clear();
    sink_set_idx.clear();
}

void SinkSetSearch::search() {
    timer.reset();

    kickstart();
    int counter = 0;

    while (!pq.empty()) {
        const size_t t(pq.top());
        dequeue();
        counter++;
        expand(vsir[t].variables);
    }

//    cerr << __PRETTY_FUNCTION__ << endl;
//    cerr << "# enqueued sets: " << counter << endl;
//    cerr << "# registered sets: " << vsir.size() << endl;
}

void SinkSetSearch::kickstart() {
    /* put singleton sets into the priority queue */
    for (int i = 0; i < static_cast<int>(task_proxy.get_variables().size()); i++) {
        enqueue(mst::singleton(i));
    }

    /* put mutex pair clique sets into the priority queue */
    mutex_pair_vars = mst::get_mutex_pairs_var();

    mutex_pair_relation = mst::get_mutex_pairs_relation();

    for (set<var_set_t>::iterator s = mutex_pair_vars.begin();
         s != mutex_pair_vars.end(); ++s) {
        enqueue(*s);
    }

    vector<set<var_set_t> > mutex_pair_cliques;
    mutex_pair_cliques.push_back(mutex_pair_vars);
    for (size_t i = 3; i <= (size_t)clique_limit; i++) {
        set<var_set_t> larger;
        if (!get_larger_mutex_pair_cliques(
                mutex_pair_cliques.back(), larger)) {
            /* if no larger cliques found, break */
            break;
        }
        mutex_pair_cliques.push_back(larger);

        for (set<var_set_t>::iterator s = larger.begin();
             s != larger.end(); ++s) {
            enqueue(*s);
        }
    }

    /* put CG strongly connected componets into the priority queue */
//    cerr << "strongly connected components :" << endl;
    vector<vector<int> > cg_succ;
    for (int i = 0; i < static_cast<int>(task_proxy.get_variables().size()); i++) {
        cg_succ.push_back(causal_graph.get_successors(i));
    }
    vector<vector<int> > sccs = SCC(cg_succ).get_result();

    for (size_t i = 0; i < sccs.size(); ++i) {
        var_set_t scc(sccs[i].begin(), sccs[i].end());
        enqueue(scc);
    }
}

bool SinkSetSearch::prune(const var_set_t &, const var_set_t &) {
    if (opt_prune == EnumPrune::CGWC_MUTEX) {
        assert(false);
//        if (is_cg_weakly_connected(S, A))
//            return false;
//        if (miasm_abstraction->mutex_control->get_bound() &&
//            is_mutex_pair_connected(S, A))
//            return false;
    }
    return true;
}

void SinkSetSearch::update_gain(const var_set_t &S, const size_t Si) {
    if (vsir[Si].ratio >= 1)
        return;

    vector<var_set_t> check_sets;

    if (opt_gain == EnumGain::ALL_ACCUR || opt_gain == EnumGain::ALL_GUESS) {
        for (size_t k = 1; k <= S.size() / 2; ++k) {
            vector<var_set_t> check_sets_k;
            k_subsets(S, k, check_sets_k);
            check_sets.insert(check_sets.end(),
                              check_sets_k.begin(),
                              check_sets_k.end());
        }
    }

    double min_ratio = 1;
//        VarSet min_K[2];
//        double min_K_r[2];
//        cerr << S << " = " << endl;

//    cerr << "gain update for" << S << endl;
    for (size_t i = 0; i < check_sets.size(); ++i) {
        var_set_t K[2];
        double K_r[2];
        K[0] = check_sets[i];
        std::set_difference(
            S.begin(), S.end(),
            K[0].begin(), K[0].end(),
            std::inserter(K[1], K[1].begin()));
//            cerr << K[0] << " * " << K[1] << endl;
        for (size_t j = 0; j < 2; j++) {
            K_r[j] = 1;
            if (vsir.contain(K[j])) {
                K_r[j] = vsir[K[j]].ratio;
            } else if (opt_gain == EnumGain::POOL_ACCUR ||
                       opt_gain == EnumGain::ALL_ACCUR) {

                compute_varset_info(K[j]);
                assert(vsir.contain(K[j]));
                K_r[j] = vsir[K[j]].ratio;
            }
        }

        if (K_r[0] * K_r[1] < min_ratio) {
            min_ratio = K_r[0] * K_r[1];
//                min_K[0] = K[0];
//                min_K[1] = K[1];
//                min_K_r[0] = K_r[0];
//                min_K_r[1] = K_r[1];
        }

        if (abs(min_ratio - vsir[Si].ratio) <
            numeric_limits<double>::epsilon()) {
            break;
        }
    }

    if (min_ratio - vsir[Si].ratio < vsir[Si].gain) {
        vsir[Si].gain = min_ratio - vsir[Si].ratio;
    }
//        cerr << min_ratio << endl
//             << subset_info[Si].ratio << endl
//             << subset_info[Si].gain << endl;

//        cerr << S << " = " << min_K[0] << " x " << min_K[1] << endl;
//        cerr << subset_info[Si].ratio
//             << " = " << min_K_r[0] << " x " << min_K_r[1] << endl;
}

bool SinkSetSearch::enqueue(const var_set_t &S, pair<size_t, size_t> P) {
    if (S.size() != 1 && time_limit_exceeded()) {
        return false;
    }

    /* duplicate check */
    if (visited.count(S)) {
        return false;
    } else {
        visited.insert(S);
    }

    size_t estimated_size = 0;

    if (P.first != numeric_limits<size_t>::max() &&
        P.second != numeric_limits<size_t>::max()) {
        const VarSetInfo &L = vsir[P.first];
        const VarSetInfo &R = vsir[P.second];
        estimated_size += combinatorial_size(L.variables, task_proxy) *
                          L.ratio;
        estimated_size += combinatorial_size(R.variables, task_proxy) *
                          R.ratio;
    } else {
        estimated_size = combinatorial_size(S, task_proxy);
    }

    /* if the abstraction is too large, then skip it */
    if (estimated_size > (size_t)size_limit) {
        return false;
    }

    compute_varset_info(S, P);

    assert(vsir.contain(S));
    const size_t Si = vsir.idx(S);

    update_gain(S, Si);

    /* do the push after the abstraction has been checked */
    pq.push(Si);

    /* check if the subset satisfies the properties being a sink set */
    if (S.size() == 1 || vsir[Si].gain > 0) {
        sink_set_idx.push_back(Si);
    }

    return true;
}

void SinkSetSearch::dequeue() {
    pq.pop();
}

void SinkSetSearch::expand(const var_set_t S) {
    if (opt_expa == EnumExpand::SINGLE) {
        for (int v = 0; v < static_cast<int>(task_proxy.get_variables().size()); ++v) {
            /* must be a new variable */
            if (S.count(v))
                continue;

            var_set_t G(S);
            G.insert(v);
            /* duplicate check */
            if (visited.count(G))
                continue;

            var_set_t A = mst::singleton(v);

            if (prune(S, A)) {
                if (!vsir.contain(G)) {
                    size_t Gi = vsir.add(G);
                    assert(vsir.contain(S) && vsir.contain(A));
                    vsir[Gi].ratio =
                        vsir[S].ratio * vsir[A].ratio;
                    vsir[Gi].gain = 0;
                }
                /* visited mark */
                visited.insert(G);
                continue;
            } else {
                enqueue(G, make_pair<size_t, size_t>(vsir.idx(S),
                                                     vsir.idx(A)));
            }
        }
    }
}

bool SinkSetSearch::get_larger_mutex_pair_cliques(
    const set<var_set_t> &current, set<var_set_t> &larger) {
    for (set<var_set_t>::iterator ia = current.begin();
         ia != current.end(); ++ia) {
        const var_set_t &old_set = *ia;
        for (int i = 0; i < static_cast<int>(task_proxy.get_variables().size()); i++) {
            if (old_set.count(i)) {
                continue;
            }

            var_set_t new_set(old_set);
            new_set.insert(i);
            if (larger.count(new_set)) {
                continue;
            }
            size_t d = 0;
            for (var_set_t::iterator j = old_set.begin();
                 j != old_set.end(); ++j) {
                if (mutex_pair_relation[i][*j])
                    d++;
            }
            if (d == old_set.size()) {
                larger.insert(new_set);
            }
        }
    }
    return larger.size();
}

bool SinkSetSearch::has_cg_predecessor(
    const int v, const var_set_t &a) const {
    vector<int> pred = causal_graph.get_predecessors(v);

    for (size_t i = 0; i < pred.size(); ++i) {
        if (a.count(i))
            return true;
    }
    return false;
}

bool SinkSetSearch::is_cg_weakly_connected(
    const var_set_t &a, const var_set_t &b) const {
    for (var_set_t::iterator u = a.begin();
         u != a.end(); ++u) {
        if (has_cg_predecessor(*u, b)) {
            return true;
        }
    }
    for (var_set_t::iterator v = b.begin();
         v != b.end(); ++v) {
        if (has_cg_predecessor(*v, a)) {
            return true;
        }
    }
    return false;
}

bool SinkSetSearch::is_mutex_pair_connected(
    const var_set_t &A, const var_set_t &B) const {
    for (var_set_t::iterator i = A.begin(); i != A.end(); ++i) {
        for (var_set_t::iterator j = B.begin(); j != B.end(); ++j) {
            var_set_t T;
            T.insert(*i);
            T.insert(*j);
            if (mutex_pair_vars.count(T))
                return true;
        }
    }
    return false;
}


void SinkSetSearch::compute_varset_info(const var_set_t &S,
                                        pair<size_t, size_t> P) {
    if (vsir.contain(S))
        return;

    VarSetInfo &vsi = vsir[vsir.add(S)];

    vsi.parent = P;

//    cerr << "compute the abstraction on " << S << endl;

    vector<var_set_t> newly_built;
    int ts_index =
        miasm_abstraction->build_transition_system(S, newly_built, vsir);
    const TransitionSystem &ts = miasm_abstraction->fts->get_ts(ts_index);

    /* initialize the ratio and gain */
    const vector<int> &var_id_set = ts.get_incorporated_variables();
    set<int> var_id_set2(var_id_set.begin(), var_id_set.end());
    assert(S == var_id_set2);
    vsi.ratio = (double)ts.get_size() /
                combinatorial_size(var_id_set2, task_proxy);
    /* defaul gain */
    vsi.gain = 1.0 - vsi.ratio;

    if (P.first != numeric_limits<size_t>::max() &&
        P.second != numeric_limits<size_t>::max()) {
        const VarSetInfo &L = vsir[P.first];
        const VarSetInfo &R = vsir[P.second];
        vsi.gain = L.ratio * R.ratio - vsi.ratio;
    }

    /* todo, the getCurrent get the RSS "recently" used, not the
     * real current usage of memory */
    if (memory_limit_exceeded()) {
        for (size_t i = 0; i < newly_built.size(); ++i) {
            if (newly_built[i].size() == 1) {
                // Do not delete cached atomic abstractions, as we cannot
                // recompute single atomic abstractions.
                continue;
            }
            miasm_abstraction->release_cache(newly_built[i]);
        }
    }
}

void SinkSetSearch::k_subsets(const var_set_t &S, const size_t k,
                              vector<var_set_t> &ks) {
    vector<vector<int> > ac = generate_combinations(S.size(), k);
    vector<int> s2v(S.begin(), S.end());
    for (size_t i = 0; i < ac.size(); i++) {
        var_set_t subset;
        for (size_t j = 0; j < k; j++) {
            subset.insert(s2v[ac[i][j]]);
        }
        /* TODO deal with complement overlapping? */
        ks.push_back(subset);
    }
}


vector<vector<int> > SinkSetSearch::generate_combinations(
    const int n, const int t) {
    assert(t > 0);
    vector<vector<int> > ac;
    vector<int> c;
    for (int i = 0; i < t; i++) {
        c.push_back(i);
    }
    int x, j;
    j = t - 1;
    ac.push_back(c);
    x = j + 1;
    while (j < t) {
        if (j >= 0) {
            c[j] = x;
            j--;
            ac.push_back(c);
            x = j + 1;
        } else {
            if ((c.size() == 1 && c[0] + 1 < n) ||
                (c.size() > 1 && c[0] + 1 < c[1])) {
                c[0]++;
                ac.push_back(c);
            } else {
                do {
                    j++;
                    c[j] = j;
                    if (j + 1 < t) {
                        x = c[j + 1] + 1;
                    }
                } while ((j + 1 == t - 1 && x == n)
                         || (j + 1 < t - 1 && x == c[j + 2]));
                j++;
            }
        }
    }
    /* TODO if we want to store this, clear unused reserved memory */
    return ac;
}

const VarSetInfoRegistry *SinkSetSearch::get_vsir() {
    return &vsir;
}

ComparatorSTLPriorityQueue::ComparatorSTLPriorityQueue(
    const shared_ptr<AbstractTask> task,
    const VarSetInfoRegistry *vsir_, const EnumPriority *priority_)
    : ComparatorVarSet(task, vsir_),
      priority(priority_) {
    /* TODO: can we specify cmp_type.e in the initilization list
     * after knowing p_order? */
    if (priority) {
        if (*priority == EnumPriority::FIFO) {
            cmp_type.e = VarSetCmpType::BY_INDEX;
        } else if (*priority == EnumPriority::RATIO) {
            cmp_type.e = VarSetCmpType::BY_RATIO;
        }
    }
}

ComparatorSTLPriorityQueue::~ComparatorSTLPriorityQueue() {
}

bool ComparatorSTLPriorityQueue::operator()(const size_t i,
                                            const size_t j) const {
    assert(vsir && priority);
    return !(ComparatorVarSet::operator ()(i, j));
}
}
