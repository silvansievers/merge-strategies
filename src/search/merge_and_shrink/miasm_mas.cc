#include "miasm_mas.h"

#include "labels.h"
//#include "merge_and_shrink_heuristic.h"
#include "merge_strategy.h"
#include "shrink_strategy.h"
#include "subset_info.h"
#include "transition_system.h"

#include "../option_parser.h"
#include "../plugin.h"
#include "../utilities.h"

#include <cassert>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

using namespace std;
using namespace mst;

MiasmAbstraction::MiasmAbstraction(const Options &opts)
    : task(get_task_from_options(opts)),
      task_proxy(*task),
      merge_strategy(opts.get<MergeStrategy *>("merge_strategy")),
      shrink_strategy(opts.get<ShrinkStrategy *>("shrink_strategy")),
      labels(opts.get<Labels *>("label_reduction")) {
    merge_strategy->initialize(task);
    labels->initialize(task_proxy);
}

MiasmAbstraction::~MiasmAbstraction() {
    delete merge_strategy;
    delete shrink_strategy;
    delete labels;
}

string MiasmAbstraction::option_key() {
    return "abstraction";
}

string MiasmAbstraction::plugin_key() {
    return "merge_and_shrink";
}

void MiasmAbstraction::release_cache(const var_set_t &var_set) {
    cerr << __PRETTY_FUNCTION__ << endl;
    assert(cache.count(var_set));
    delete cache[var_set];
    cache.erase(var_set);
}

void MiasmAbstraction::release_cache() {
    cerr << __PRETTY_FUNCTION__ << endl;
    for (map<var_set_t, TransitionSystem *>::iterator i = cache.begin();
         i != cache.end(); ++i) {
//        cerr << i->first << endl;
        delete i->second;
    }
    map<var_set_t, TransitionSystem *>().swap(cache);
}

TransitionSystem *MiasmAbstraction::build_transition_system(
    const var_set_t &G, vector<var_set_t> &newly_built,
    const VarSetInfoRegistry &vsir) {
    if (cache.count(G)) {
        cerr << "old: " << G << endl;
        return cache[G];
    }
    assert(!G.empty());
    /* will do once only */
    if (G.size() == 1) {
        vector<TransitionSystem *> atomic;
        TransitionSystem::build_atomic_transition_systems(task_proxy, atomic, labels);

        /* remove the atomic abstraction if its variable is not involved */
        for (var_t i = 0; (size_t)i < atomic.size(); ++i) {
            var_set_t s = mst::singleton(i);
            assert(!cache.count(s));
            newly_built.push_back(s);
            cerr << "new: " << s << endl;
            cache.insert(pair<var_set_t, TransitionSystem *>(s, atomic[i]));
        }

        assert(cache.count(G));
        return cache[G];
    }

    var_set_t left_set, right_set;

    if (vsir.contain(G)) {
        size_t pl = vsir[G].parent.first;
        size_t pr = vsir[G].parent.second;
        if (pl != numeric_limits<size_t>::max() &&
            pr != numeric_limits<size_t>::max()) {
            left_set = vsir[pl].variables;
            right_set = vsir[pr].variables;
        }
    }

    if (left_set.empty() && right_set.empty()) {
        // TODO: currently the abstraction on varset is constructed by simply
        // merging corresponding atomic abstractions in the default order
        // without any shrinking or label reduction
        // re-implement this to allow arbitrary merging, shrinking and label reduction
        vector<var_t> ordered(G.begin(), G.end());
        left_set = set<var_t>(ordered.begin(), ordered.end() - 1);
        right_set.insert(ordered.back());

//        if (!vsir.contain(left_set))
//            vsir.add(left_set);
//        if (!vsir.contain(right_set))
//            vsir.add(right_set);

//        if (!vsir.contain(G)) {
//            vsir.add(G);
//            vsir[G].parent = make_pair<size_t, size_t>(
//                vsir.idx(left_set), vsir.idx(left_set));
//        }
    }

    cerr << left_set << ", " << right_set << endl;


    TransitionSystem *left = build_transition_system(left_set,
                                                     newly_built, vsir);
    TransitionSystem *right = build_transition_system(right_set,
                                                      newly_built, vsir);

    TransitionSystem *root = new CompositeTransitionSystem(task_proxy, labels, left, right);

    newly_built.push_back(G);
    cache.insert(pair<var_set_t, TransitionSystem *>(G, root));
    assert(cache.count(G));
    cerr << "new: " << G << endl;
    return cache[G];
}

static MiasmAbstraction *_parse(OptionParser &parser) {
    // Merge strategy option.
    parser.add_option<MergeStrategy *>(
        "merge_strategy",
        "merge strategy; choose between merge_linear with various variable "
        "orderings and merge_dfp.");

    // Shrink strategy option.
    parser.add_option<ShrinkStrategy *>(
        "shrink_strategy",
        "shrink strategy; choose between shrink_fh and shrink_bisimulation. "
        "A good configuration for bisimulation based shrinking is: "
        "shrink_bisimulation(max_states=50000, max_states_before_merge=50000, "
        "threshold=1, greedy=false)");

    // Label reduction option.
    parser.add_option<Labels *>("label_reduction",
                                "Choose relevant options for label reduction. "
                                "Also note the interaction with shrink strategies.");

    // For AbstractTask
    Heuristic::add_options_to_parser(parser);
    Options opts = parser.parse();

    if (parser.dry_run()) {
        return 0;
    } else {
        return new MiasmAbstraction(opts);
    }
}

static Plugin<MiasmAbstraction> _plugin(MiasmAbstraction::plugin_key(), _parse);

