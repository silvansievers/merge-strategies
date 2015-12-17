#include "label_equivalence_relation.h"

#include "labels.h"

//#include "../utilities.h"

#include <cassert>
#include <iostream>


using namespace std;

namespace MergeAndShrink {
LabelEquivalenceRelation::LabelEquivalenceRelation(const Labels &labels)
    : labels(labels) {
    grouped_labels.reserve(labels.get_max_size());
    label_to_positions.resize(labels.get_max_size());
}

LabelEquivalenceRelation::LabelEquivalenceRelation(const LabelEquivalenceRelation &other)
    : labels(other.labels),
      label_to_positions(other.label_to_positions) {
    // we need to reserve max size to ensure that no move occurs in grouped_labels
    // otherwise, iterators to elements of list<int> of LabelGroup could become
    // invalid!
    grouped_labels.reserve(labels.get_max_size());
    for (size_t group_id = 0; group_id < other.grouped_labels.size(); ++group_id) {
        const LabelGroup &label_group = other.grouped_labels[group_id];
        grouped_labels.push_back(LabelGroup(label_group));
        // TODO: change so that we do not first create the group and then update
        // the iterators.
        LabelGroup &new_group = grouped_labels.back();
        // we also need to update label_to_positions with correct iterators
        for (LabelIter label_it = new_group.begin();
             label_it != new_group.end(); ++label_it) {
            label_to_positions[*label_it] = make_pair(group_id, label_it);
        }
    }
    assert(grouped_labels == other.grouped_labels);
    assert(*this == other);
}

void LabelEquivalenceRelation::add_label_to_group(int group_id,
                                                  int label_no) {
    LabelIter label_it = grouped_labels[group_id].insert(label_no);
    assert(*label_it == label_no);
    label_to_positions[label_no] = make_pair(group_id, label_it);

    int label_cost = labels.get_label_cost(label_no);
    if (label_cost < grouped_labels[group_id].get_cost())
        grouped_labels[group_id].set_cost(label_cost);
}

void LabelEquivalenceRelation::apply_label_mapping(
    const vector<pair<int, vector<int>>> &label_mapping,
    const unordered_set<int> *affected_group_ids) {
    for (const pair<int, vector<int>> &mapping : label_mapping) {
        int new_label_no = mapping.first;
        const vector<int> &old_label_nos = mapping.second;

        // Add new label to group
        int canonical_group_id = get_group_id(old_label_nos.front());
        if (!affected_group_ids) {
            add_label_to_group(canonical_group_id, new_label_no);
        } else {
            add_label_group({new_label_no});
        }

        // Remove old labels from group
        for (int old_label_no : old_label_nos) {
            if (!affected_group_ids) {
                assert(canonical_group_id == get_group_id(old_label_no));
            }
            LabelIter label_it = label_to_positions[old_label_no].second;
            grouped_labels[canonical_group_id].erase(label_it);
        }
    }

    if (affected_group_ids) {
        // Recompute the cost of all affected label groups.
        const unordered_set<int> &group_ids = *affected_group_ids;
        for (int group_id : group_ids) {
            LabelGroup &label_group = grouped_labels[group_id];
            // Setting cost to infinity for empty groups does not hurt.
            label_group.set_cost(INF);
            for (int label_no : label_group) {
                int cost = labels.get_label_cost(label_no);
                if (cost < label_group.get_cost()) {
                    label_group.set_cost(cost);
                }
            }
        }
    }
}

void LabelEquivalenceRelation::move_group_into_group(
    int from_group_id, int to_group_id) {
    assert(!is_empty_group(from_group_id));
    assert(!is_empty_group(to_group_id));
    LabelGroup &from_group = grouped_labels[from_group_id];
    for (int label_no : from_group) {
        add_label_to_group(to_group_id, label_no);
    }
    from_group.clear();
}

int LabelEquivalenceRelation::add_label_group(const vector<int> &new_labels) {
    int new_group_id = grouped_labels.size();
    grouped_labels.push_back(LabelGroup());
    for (int label_no : new_labels) {
        add_label_to_group(new_group_id, label_no);
    }
    return new_group_id;
}

bool LabelEquivalenceRelation::consistent() {
    bool consistent = true;
    for (size_t group_id = 0; group_id < grouped_labels.size(); ++group_id) {
        LabelGroup &label_group = grouped_labels[group_id];
        for (LabelIter label_it = label_group.begin();
             label_it != label_group.end(); ++label_it) {
            int label_no = *label_it;
            if (label_to_positions[label_no].first != static_cast<int>(group_id) ||
                *label_to_positions[label_no].second != *label_it) {
                cout << label_no << " has inconsistent entries" << endl;
                cout << label_to_positions[label_no].first << " "
                     << group_id << endl;
                cout << *label_to_positions[label_no].second << " "
                     << *label_it << endl;
                consistent = false;
                break;
            }
        }
    }
    return consistent;
}

bool LabelEquivalenceRelation::operator==(const LabelEquivalenceRelation &other) const {
    assert(grouped_labels == other.grouped_labels);
    assert(label_to_positions.size() == other.label_to_positions.size());
    bool label_to_positions_consistent = true;
    for (size_t i = 0; i < label_to_positions.size(); ++i) {
        if (label_to_positions[i].first != other.label_to_positions[i].first) {
            // we cannot compare the iterators, i.e. the 'second' entries,
            // because they may be invalid
            label_to_positions_consistent = false;
            break;
        }
    }
    assert(label_to_positions_consistent);
    return grouped_labels == other.grouped_labels && label_to_positions_consistent;
}

void LabelEquivalenceRelation::dump() const {
    cout << "label equiv rel:" << endl;
    for (const LabelGroup &group : grouped_labels) {
        if (!group.empty()) {
            for (LabelConstIter label_it = group.begin(); label_it != group.end(); ++label_it) {
                cout << *label_it << ", ";
            }
            cout << endl;
        }
    }
}
}
