#include "shrink_strategy.h"

#include "factored_transition_system.h"
#include "transition_system.h"

#include "../plugin.h"

#include <iostream>

using namespace std;

namespace merge_and_shrink {
ShrinkStrategy::ShrinkStrategy() {
}

ShrinkStrategy::~ShrinkStrategy() {
}

bool ShrinkStrategy::shrink(
    FactoredTransitionSystem &fts,
    int index,
    int target,
    bool silent) {
    StateEquivalenceRelation equivalence_relation;
    compute_equivalence_relation(fts, index, target, equivalence_relation);
    // TODO: We currently violate this; see issue250
    //assert(equivalence_relation.size() <= new_size);
    return fts.apply_abstraction(index, equivalence_relation, silent);
}

int ShrinkStrategy::compute_size_after_perfect_shrink(
    const FactoredTransitionSystem &fts,
    int index) {
    StateEquivalenceRelation equivalence_relation;
    compute_equivalence_relation(fts, index, fts.get_ts(index).get_size(), equivalence_relation);
    return equivalence_relation.size();
}

void ShrinkStrategy::dump_options() const {
    cout << "Shrink strategy options: " << endl;
    cout << "Type: " << name() << endl;
    dump_strategy_specific_options();
}

string ShrinkStrategy::get_name() const {
    return name();
}

static PluginTypePlugin<ShrinkStrategy> _type_plugin(
    "ShrinkStrategy",
    "This page describes the various shrink strategies supported "
    "by the planner.");
}
