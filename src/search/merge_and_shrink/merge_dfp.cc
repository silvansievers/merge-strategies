#include "merge_dfp.h"

#include "labels.h"
#include "transition_system.h"

#include "../option_parser.h"
#include "../plugin.h"

#include <algorithm>
#include <cassert>
#include <iostream>

using namespace std;


MergeDFP::MergeDFP(const Options &options)
    : MergeStrategy(), use_cost(options.get<bool>("use_cost")) {
}

void MergeDFP::initialize(const shared_ptr<AbstractTask> task) {
    MergeStrategy::initialize(task);
    /*
      n := remaining_merges + 1 is the number of variables of the planning task
      and thus the number of atomic transition systems. These will be stored at
      indices 0 to n-1 and thus n is the index at which the first composite
      transition system will be stored at.
    */
    border_atomics_composites = remaining_merges + 1;
}

int MergeDFP::get_corrected_index(int index) const {
    /*
      This method assumes that we iterate over the vector of all
      transition systems in inverted order (from back to front). It returns the
      unmodified index as long as we are in the range of composite
      transition systems (these are thus traversed in order from the last one
      to the first one) and modifies the index otherwise so that the order
      in which atomic transition systems are considered is from the first to the
      last one (from front to back). This is to emulate the previous behavior
      when new transition systems were not inserted after existing transition systems,
      but rather replaced arbitrarily one of the two original transition systems.
    */
    assert(index >= 0);
    if (index >= border_atomics_composites)
        return index;
    return border_atomics_composites - 1 - index;
}

void MergeDFP::compute_label_ranks(const TransitionSystem *transition_system,
                                   vector<int> &label_ranks) const {
    int num_labels = transition_system->get_num_labels();
    // Irrelevant (and inactive, i.e. reduced) labels have a dummy rank of -1
    label_ranks.resize(num_labels, -1);

    for (TSConstIterator group_it = transition_system->begin();
         group_it != transition_system->end(); ++group_it) {
        // Relevant labels with no transitions have a rank of infinity.
        int label_rank = INF;
        const vector<Transition> &transitions = group_it.get_transitions();
        bool group_relevant = false;
        if (static_cast<int>(transitions.size()) == transition_system->get_size()) {
            /*
              A label group is irrelevant in the earlier notion if it has
              exactly a self loop transition for every state.
            */
            for (size_t i = 0; i < transitions.size(); ++i) {
                if (transitions[i].target != transitions[i].src) {
                    group_relevant = true;
                    break;
                }
            }
        } else {
            group_relevant = true;
        }
        if (!group_relevant) {
            label_rank = -1;
        } else {
            for (size_t i = 0; i < transitions.size(); ++i) {
                const Transition &t = transitions[i];
                label_rank = min(label_rank, transition_system->get_goal_distance(t.target));
            }
        }
        for (LabelConstIter label_it = group_it.begin();
             label_it != group_it.end(); ++label_it) {
            int label_no = *label_it;
            if (label_rank < label_ranks[label_no] || label_ranks[label_no] == -1) {
                label_ranks[label_no] = label_rank;
            }
        }
    }
}

void MergeDFP::compute_relevant_labels(const TransitionSystem *ts,
                                       vector<bool> &relevant_labels) const {
    // relevant_labels[i] is true for label i iff it is relevant and not dead
    // in transition system ts.
    assert(relevant_labels.empty());
    int num_labels = ts->get_num_labels();
    relevant_labels.resize(num_labels, true);
    for (TSConstIterator group_it = ts->begin();
         group_it != ts->end(); ++group_it) {
        const vector<Transition> &transitions = group_it.get_transitions();
        bool group_relevant = false;
        if (static_cast<int>(transitions.size()) == ts->get_size()) {
            /*
              A label group is irrelevant in the earlier notion if it has
              exactly a self loop transition for every state.
            */
            for (size_t i = 0; i < transitions.size(); ++i) {
                if (transitions[i].target != transitions[i].src) {
                    group_relevant = true;
                    break;
                }
            }
        } else if (!transitions.empty()){
            group_relevant = true;
        }
        if (!group_relevant) {
            for (LabelConstIter label_it = group_it.begin();
                 label_it != group_it.end(); ++label_it) {
                relevant_labels[*label_it] = false;
            }
        }
    }
}

pair<int, int> MergeDFP::get_next_cost(const vector<TransitionSystem *> &all_transition_systems) {
    TransitionSystem *some_ts = 0;
    for (TransitionSystem *ts : all_transition_systems) {
        if (ts) {
            some_ts = ts;
            break;
        }
    }
    vector<vector<bool> > transition_system_relevant_labels;
    transition_system_relevant_labels.resize(all_transition_systems.size());
    const shared_ptr<Labels> labels = some_ts->get_labels();
    int next_index1 = -1;
    int next_index2 = -1;
    int max_cost = -1;
    for (size_t i = 0; i < all_transition_systems.size(); ++i) {
        TransitionSystem *ts1 = all_transition_systems[i];
        if (ts1) {
            vector<bool> &relevant_labels1 = transition_system_relevant_labels[i];
            if (relevant_labels1.empty()) {
                compute_relevant_labels(ts1, relevant_labels1);
            }
            for (size_t j = i + 1; j < all_transition_systems.size(); ++j) {
                TransitionSystem *ts2 = all_transition_systems[j];
                if (ts2) {
                    vector<bool> &relevant_labels2 = transition_system_relevant_labels[j];
                    if (relevant_labels2.empty()) {
                        compute_relevant_labels(ts2, relevant_labels2);
                    }
                    int num_labels = labels->get_size();
                    int label_costs = 0;
                    for (int label_no = 0; label_no < num_labels; ++label_no) {
                        if (labels->is_current_label(label_no)) {
                            if (relevant_labels1[label_no] || relevant_labels2[label_no]) {
                                label_costs += labels->get_label_cost(label_no);
                            }
                        }
                    }
//                    cout << ts1->tag() << ts2->tag() << "cost: " << label_costs << endl;
                    if (label_costs > max_cost) {
                        max_cost = label_costs;
                        next_index1 = i;
                        next_index2 = j;
                    }
                }
            }
        }
    }
    assert(next_index1 != -1);
    assert(next_index2 != -1);
    cout << "Next pair of indices: (" << next_index1 << ", " << next_index2 << ")" << endl;
    --remaining_merges;
    return make_pair(next_index1, next_index2);
}

pair<int, int> MergeDFP::get_next(const vector<TransitionSystem *> &all_transition_systems) {
    assert(initialized());
    assert(!done());

    if (use_cost) {
        return get_next_cost(all_transition_systems);
    }

    vector<const TransitionSystem *> sorted_transition_systems;
    vector<int> indices_mapping;
    vector<vector<int> > transition_system_label_ranks;
    /*
      Precompute a vector sorted_transition_systems which contains all exisiting
      transition systems from all_transition_systems in the desired order and
      compute label ranks.
    */
    for (int i = all_transition_systems.size() - 1; i >= 0; --i) {
        /*
          We iterate from back to front, considering the composite
          transition systems in the order from "most recently added" (= at the back
          of the vector) to "first added" (= at border_atomics_composites).
          Afterwards, we consider the atomic transition systems in the "regular"
          order from the first one until the last one. See also explanation
          at get_corrected_index().
        */
        int ts_index = get_corrected_index(i);
        const TransitionSystem *transition_system = all_transition_systems[ts_index];
        if (transition_system) {
            sorted_transition_systems.push_back(transition_system);
            indices_mapping.push_back(ts_index);
            transition_system_label_ranks.push_back(vector<int>());
            vector<int> &label_ranks = transition_system_label_ranks.back();
            compute_label_ranks(transition_system, label_ranks);
        }
    }

    int next_index1 = -1;
    int next_index2 = -1;
    int minimum_weight = INF;
    for (size_t i = 0; i < sorted_transition_systems.size(); ++i) {
        const TransitionSystem *transition_system1 = sorted_transition_systems[i];
        assert(transition_system1);
        const vector<int> &label_ranks1 = transition_system_label_ranks[i];
        assert(!label_ranks1.empty());
        for (size_t j = i + 1; j < sorted_transition_systems.size(); ++j) {
            const TransitionSystem *transition_system2 = sorted_transition_systems[j];
            assert(transition_system2);
            if (transition_system1->is_goal_relevant()
                || transition_system2->is_goal_relevant()) {
                vector<int> &label_ranks2 = transition_system_label_ranks[j];
                assert(!label_ranks2.empty());
                assert(label_ranks1.size() == label_ranks2.size());
                int pair_weight = INF;
                for (size_t k = 0; k < label_ranks1.size(); ++k) {
                    if (label_ranks1[k] != -1 && label_ranks2[k] != -1) {
                        // label is relevant in both transition_systems
                        int max_label_rank = max(label_ranks1[k], label_ranks2[k]);
                        pair_weight = min(pair_weight, max_label_rank);
                    }
                }
                if (pair_weight < minimum_weight) {
                    minimum_weight = pair_weight;
                    next_index1 = indices_mapping[i];
                    next_index2 = indices_mapping[j];
                    assert(all_transition_systems[next_index1] == transition_system1);
                    assert(all_transition_systems[next_index2] == transition_system2);
                }
            }
        }
    }
    if (next_index1 == -1) {
        /*
          No pair with finite weight has been found. In this case, we simply
          take the first pair according to our ordering consisting of at
          least one goal relevant transition system.
        */
        assert(next_index2 == -1);
        assert(minimum_weight == INF);

        for (size_t i = 0; i < sorted_transition_systems.size(); ++i) {
            const TransitionSystem *transition_system1 = sorted_transition_systems[i];
            for (size_t j = i + 1; j < sorted_transition_systems.size(); ++j) {
                const TransitionSystem *transition_system2 = sorted_transition_systems[j];
                if (transition_system1->is_goal_relevant()
                    || transition_system2->is_goal_relevant()) {
                    next_index1 = indices_mapping[i];
                    next_index2 = indices_mapping[j];
                    assert(all_transition_systems[next_index1] == transition_system1);
                    assert(all_transition_systems[next_index2] == transition_system2);
                }
            }
        }
    }
    /*
      There always exists at least one goal relevant transition system,
      assuming that the global goal specification is non-empty. Hence at
      this point, we must have found a pair of transition systems to merge.
    */
    // NOT true if used on a subset of transitions!
    if (next_index1 == -1 || next_index2 == -1) {
        assert(next_index1 == -1 && next_index2 == -1);
        assert(minimum_weight == INF);
        size_t ts_index = 0;
        size_t other_ts_index = 1;
        assert(sorted_transition_systems[ts_index]);
        assert(sorted_transition_systems[other_ts_index]);
        next_index1 = indices_mapping[ts_index];
        next_index2 = indices_mapping[other_ts_index];
        assert(all_transition_systems[next_index1] == sorted_transition_systems[ts_index]);
        assert(all_transition_systems[next_index2] == sorted_transition_systems[other_ts_index]);
    }
    assert(next_index1 != -1);
    assert(next_index2 != -1);
    cout << "Next pair of indices: (" << next_index1 << ", " << next_index2 << ")" << endl;
    --remaining_merges;
    return make_pair(next_index1, next_index2);
}

string MergeDFP::name() const {
    return "dfp";
}

static shared_ptr<MergeStrategy>_parse(OptionParser &parser) {
    parser.document_synopsis(
        "Merge strategy DFP",
        "This merge strategy implements the algorithm originally described in the "
        "paper \"Directed model checking with distance-preserving abstractions\" "
        "by Draeger, Finkbeiner and Podelski (SPIN 2006), adapted to planning in "
        "the following paper:\n\n"
        " * Silvan Sievers, Martin Wehrle, and Malte Helmert.<<BR>>\n"
        " [Generalized Label Reduction for Merge-and-Shrink Heuristics "
        "http://ai.cs.unibas.ch/papers/sievers-et-al-aaai2014.pdf].<<BR>>\n "
        "In //Proceedings of the 28th AAAI Conference on Artificial "
        "Intelligence (AAAI 2014)//, pp. 2358-2366. AAAI Press 2014.");
    parser.add_option<bool>("use_cost", "use a different merge strategy", "false");
    Options options = parser.parse();
    if (parser.dry_run())
        return nullptr;
    else
        return make_shared<MergeDFP>(options);
}

static PluginShared<MergeStrategy> _plugin("merge_dfp", _parse);
