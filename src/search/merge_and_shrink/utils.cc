#include "utils.h"

#include "distances.h"
#include "factored_transition_system.h"
#include "label_equivalence_relation.h"
#include "labels.h"
#include "shrink_strategy.h"
#include "transition_system.h"

#include "../utils/math.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <unordered_map>

using namespace std;

namespace merge_and_shrink {
pair<int, int> compute_shrink_sizes(
    int size1,
    int size2,
    int max_states_before_merge,
    int max_states_after_merge) {
    // Bound both sizes by max allowed size before merge.
    int new_size1 = min(size1, max_states_before_merge);
    int new_size2 = min(size2, max_states_before_merge);

    if (!utils::is_product_within_limit(
            new_size1, new_size2, max_states_after_merge)) {
        int balanced_size = int(sqrt(max_states_after_merge));

        if (new_size1 <= balanced_size) {
            // Size of the first transition system is small enough. Use whatever
            // is left for the second transition system.
            new_size2 = max_states_after_merge / new_size1;
        } else if (new_size2 <= balanced_size) {
            // Inverted case as before.
            new_size1 = max_states_after_merge / new_size2;
        } else {
            // Both transition systems are too big. We set both target sizes
            // to balanced_size. An alternative would be to set one to
            // N1 = balanced_size and the other to N2 = max_states_after_merge /
            // balanced_size, to get closer to the allowed maximum.
            // However, this would make little difference (N2 would
            // always be N1, N1 + 1 or N1 + 2), and our solution has the
            // advantage of treating the transition systems symmetrically.
            new_size1 = balanced_size;
            new_size2 = balanced_size;
        }
    }
    assert(new_size1 <= size1 && new_size2 <= size2);
    assert(new_size1 <= max_states_before_merge);
    assert(new_size2 <= max_states_before_merge);
    assert(new_size1 * new_size2 <= max_states_after_merge);
    return make_pair(new_size1, new_size2);
}

bool shrink_transition_system(
    FactoredTransitionSystem &fts,
    int index,
    int new_size,
    int shrink_threshold_before_merge,
    const ShrinkStrategy &shrink_strategy,
    Verbosity verbosity) {
    const TransitionSystem &ts = fts.get_ts(index);
    assert(ts.is_solvable());
    int num_states = ts.get_size();
    if (num_states > min(new_size, shrink_threshold_before_merge)) {
        if (verbosity >= Verbosity::VERBOSE) {
            cout << ts.tag() << "current size: " << num_states;
            if (new_size < num_states)
                cout << " (new size limit: " << new_size;
            else
                cout << " (shrink threshold: " << shrink_threshold_before_merge;
            cout << ")" << endl;
        }
        return shrink_strategy.shrink(fts, index, new_size, verbosity);
    }
    return false;
}

bool is_goal_relevant(const TransitionSystem &ts) {
    int num_states = ts.get_size();
    for (int state = 0; state < num_states; ++state) {
        if (!ts.is_goal_state(state)) {
            return true;
        }
    }
    return false;
}

int shrink_and_merge_temporarily(
    FactoredTransitionSystem &fts,
    int ts_index1,
    int ts_index2,
    const ShrinkStrategy &shrink_strategy,
    int max_states,
    int max_states_before_merge,
    int shrink_threshold_before_merge) {
    // Copy the transition systems (distances etc)
    int copy_ts_index1 = fts.copy(ts_index1);
    int copy_ts_index2 = fts.copy(ts_index2);
    pair<int, int> shrink_sizes =
        compute_shrink_sizes(fts.get_ts(copy_ts_index1).get_size(),
                             fts.get_ts(copy_ts_index2).get_size(),
                             max_states,
                             max_states_before_merge);

    // Shrink before merge
    Verbosity verbosity = Verbosity::SILENT;
    shrink_transition_system(
        fts,
        copy_ts_index1,
        shrink_sizes.first,
        shrink_threshold_before_merge,
        shrink_strategy,
        verbosity);
    shrink_transition_system(
        fts,
        copy_ts_index2,
        shrink_sizes.second,
        shrink_threshold_before_merge,
        shrink_strategy,
        verbosity);

    // Perform the merge and temporarily add it to fts
    const bool finalize_if_unsolvable = false;
    const bool invalidating_merge = true;
    int merge_index = fts.merge(
        copy_ts_index1,
        copy_ts_index2,
        verbosity,
        finalize_if_unsolvable,
        invalidating_merge
        );
    return merge_index;
}

int compute_number_of_product_transitions(
    const TransitionSystem &ts1, const TransitionSystem &ts2) {
    // TODO: this is copied from the merge constructor of TransitionSystem
    /*
      Note that this computes the number of tranistions in the product
      without considering possible shrinking due to unreachable or
      irrelevant states, which hence may reduce the actual number of
      transitions in the product.
    */
    int number_of_transitions = 0;
    for (const GroupAndTransitions &gat : ts1) {
        const LabelGroup &group1 = gat.label_group;
        const vector<Transition> &transitions1 = gat.transitions;

        // Distribute the labels of this group among the "buckets"
        // corresponding to the groups of ts2.
        unordered_map<int, vector<double>> buckets;
        for (int label_no : group1) {
            int group2_id = ts2.get_group_id_for_label(label_no);
            buckets[group2_id].push_back(label_no);
        }

        // Now create the new groups together with their transitions.
        for (const auto &bucket : buckets) {
            const vector<Transition> &transitions2 =
                ts2.get_transitions_for_group_id(bucket.first);
            int new_transitions_for_new_group = transitions1.size() * transitions2.size();
            number_of_transitions += new_transitions_for_new_group;
        }
    }
    return number_of_transitions;
}

double compute_average_h_value(const Distances &distances) {
    int num_states = distances.get_num_states();
    int sum_distances = 0;
    for (int state = 0; state < num_states; ++state) {
        sum_distances += distances.get_goal_distance(state);
    }
    if (num_states == 0) {
        // For unsolvable transition systems
        return INF;
    }
    return static_cast<double>(sum_distances) / static_cast<double>(num_states);
}

void compute_irrelevant_labels(const FactoredTransitionSystem &fts,
                               vector<vector<bool>> &ts_index_to_irrelevant_labels) {
    int num_ts = fts.get_size();
    ts_index_to_irrelevant_labels.resize(num_ts, vector<bool>());
    int num_labels = fts.get_labels().get_size();
    for (int ts_index = 0; ts_index < num_ts; ++ts_index) {
        if (fts.is_active(ts_index)) {
            vector<bool> irrelevant_labels(num_labels, false);
            const TransitionSystem &ts = fts.get_ts(ts_index);
            for (const GroupAndTransitions &gat : ts) {
                const vector<Transition> &transitions = gat.transitions;
                bool group_relevant = false;
                if (static_cast<int>(transitions.size()) == ts.get_size()) {
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
                    const LabelGroup &label_group = gat.label_group;
                    for (int label_no : label_group) {
                        irrelevant_labels[label_no] = true;
                    }
                }
            }
            ts_index_to_irrelevant_labels[ts_index] = irrelevant_labels;
        }
    }
}
}
