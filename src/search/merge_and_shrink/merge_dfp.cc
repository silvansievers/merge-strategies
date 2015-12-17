#include "merge_dfp.h"

#include "distances.h"
#include "factored_transition_system.h"
#include "label_equivalence_relation.h"
#include "transition_system.h"
#include "types.h"

#include "../globals.h"
#include "../option_parser.h"
#include "../option_parser_util.h"
#include "../plugin.h"
#include "../task_proxy.h"

#include "../utils/rng.h"

#include <algorithm>
#include <cassert>
#include <iostream>
#include <unordered_map>

using namespace std;


namespace MergeAndShrink {
MergeDFP::MergeDFP(const Options &options)
    : MergeStrategy(),
      atomic_ts_order(AtomicTSOrder(options.get_enum("atomic_ts_order"))),
      product_ts_order(ProductTSOrder(options.get_enum("product_ts_order"))),
      atomic_before_product(options.get<bool>("atomic_before_product")),
      randomized_order(options.get<bool>("randomized_order")) {
}

void MergeDFP::initialize(const shared_ptr<AbstractTask> task) {
    MergeStrategy::initialize(task);
    TaskProxy task_proxy(*task);
    compute_ts_order(task_proxy,
                     atomic_ts_order,
                     product_ts_order,
                     atomic_before_product,
                     randomized_order,
                     transition_system_order);

    if (!randomized_order && (atomic_ts_order == REGULAR &&
                              product_ts_order == NEW_TO_OLD &&
                              !atomic_before_product)) {
        int num_variables = task_proxy.get_variables().size();
        int max_transition_system_count = num_variables * 2 - 1;
        vector<int> original_dfp_order;
        for (int i = max_transition_system_count - 1; i >= 0; --i) {
            int corrected_index = i;
            if (i < num_variables) {
                corrected_index = num_variables - 1 - i;
            }
            original_dfp_order.push_back(corrected_index);
        }
        assert(transition_system_order == original_dfp_order);
    }
}

void MergeDFP::compute_label_ranks(FactoredTransitionSystem &fts,
                                   int index,
                                   vector<int> &label_ranks) const {
    const TransitionSystem &ts = fts.get_ts(index);
    const Distances &distances = fts.get_dist(index);
    int num_labels = fts.get_num_labels();
    // Irrelevant (and inactive, i.e. reduced) labels have a dummy rank of -1
    label_ranks.resize(num_labels, -1);

    for (const GroupAndTransitions &gat : ts) {
        const LabelGroup &label_group = gat.label_group;
        const vector<Transition> &transitions = gat.transitions;
        // Relevant labels with no transitions have a rank of infinity.
        int label_rank = INF;
        bool group_relevant = false;
        if (static_cast<int>(transitions.size()) == ts.get_size()) {
            /*
              A label group is irrelevant in the earlier notion if it has
              exactly a self loop transition for every state.
            */
            for (const Transition &transition : transitions) {
                if (transition.target != transition.src) {
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
            for (const Transition &transition : transitions) {
                label_rank = min(label_rank,
                                 distances.get_goal_distance(transition.target));
            }
        }
        for (int label_no : label_group) {
            label_ranks[label_no] = label_rank;
        }
    }
}

pair<int, int> MergeDFP::get_next_dfp(
    FactoredTransitionSystem &fts,
    const vector<int> &sorted_active_ts_indices) {
    int next_index1 = -1;
    int next_index2 = -1;
    int first_valid_pair_index1 = -1;
    int first_valid_pair_index2 = -1;
    int minimum_weight = INF;
    vector<vector<int>> transition_system_label_ranks(sorted_active_ts_indices.size());
    unordered_map<int, int> weight_to_count;
    // Go over all pairs of transition systems and compute their weight.
    for (size_t i = 0; i < sorted_active_ts_indices.size(); ++i) {
        int ts_index1 = sorted_active_ts_indices[i];
        assert(fts.is_active(ts_index1));
        vector<int> &label_ranks1 = transition_system_label_ranks[i];
        if (label_ranks1.empty()) {
            compute_label_ranks(fts, ts_index1, label_ranks1);
        }
        for (size_t j = i + 1; j < sorted_active_ts_indices.size(); ++j) {
            int ts_index2 = sorted_active_ts_indices[j];
            assert(fts.is_active(ts_index2));
            if (fts.get_ts(ts_index1).is_goal_relevant()
                || fts.get_ts(ts_index2).is_goal_relevant()) {
                // Only consider pairs where at least one component is goal relevant.

                // TODO: the 'old' code that took the 'first' pair in case of
                // no finite pair weight could be found, actually took the last
                // one, so we do the same here for the moment.
//                if (first_valid_pair_index1 == -1) {
                // Remember the first such pair
//                    assert(first_valid_pair_index2 == -1);
                first_valid_pair_index1 = ts_index1;
                first_valid_pair_index2 = ts_index2;
//                }

                // Compute the weight associated with this pair
                vector<int> &label_ranks2 = transition_system_label_ranks[j];
                if (label_ranks2.empty()) {
                    compute_label_ranks(fts, ts_index2, label_ranks2);
                };
                assert(label_ranks1.size() == label_ranks2.size());
                int pair_weight = INF;
                for (size_t k = 0; k < label_ranks1.size(); ++k) {
                    if (label_ranks1[k] != -1 && label_ranks2[k] != -1) {
                        // label is relevant in both transition_systems
                        int max_label_rank = max(label_ranks1[k], label_ranks2[k]);
                        pair_weight = min(pair_weight, max_label_rank);
                    }
                }
                if (!weight_to_count.count(pair_weight)) {
                    weight_to_count[pair_weight] = 1;
                } else {
                    weight_to_count[pair_weight] += 1;
                }
                if (pair_weight < minimum_weight) {
                    minimum_weight = pair_weight;
                    next_index1 = ts_index1;
                    next_index2 = ts_index2;
                }
            }
        }
    }

    if (next_index1 == -1) {
        /*
          TODO: this is not correct (see above)! we take the *last* pair.
          We should eventually change this to be a random ordering.

          No pair with finite weight has been found. In this case, we simply
          take the first pair according to our ordering consisting of at
          least one goal relevant transition system which we compute in the
          loop before. There always exists such a pair assuming that the
          global goal specification is non-empty.
        */
        assert(next_index2 == -1);
        assert(minimum_weight == INF);
        // The above text is *not* true if only called for a subset of
        // transition systems!
        if (first_valid_pair_index1 == -1) {
            assert(first_valid_pair_index2 == -1);
            next_index1 = sorted_active_ts_indices[0];
            next_index2 = sorted_active_ts_indices[1];
        } else {
            assert(first_valid_pair_index1 != -1);
            assert(first_valid_pair_index2 != -1);
            next_index1 = first_valid_pair_index1;
            next_index2 = first_valid_pair_index2;
        }
    }

    int minimum_weight_pair_count = weight_to_count[minimum_weight];
    assert(minimum_weight_pair_count >= 1);
    if (minimum_weight_pair_count > 1) {
        ++iterations_with_tiebreaking;
        total_tiebreaking_pair_count += minimum_weight_pair_count;
    }

    assert(next_index1 != -1);
    assert(next_index2 != -1);
    return make_pair(next_index1, next_index2);
}

pair<int, int> MergeDFP::get_next(FactoredTransitionSystem &fts) {
    assert(initialized());
    assert(!done());

    /*
      Precompute a vector sorted_active_ts_indices which contains all exisiting
      transition systems in the given order and compute label ranks.
    */
    assert(!transition_system_order.empty());
    vector<int> sorted_active_ts_indices;
    vector<vector<int>> transition_system_label_ranks;
    for (size_t tso_index = 0; tso_index < transition_system_order.size(); ++tso_index) {
        int ts_index = transition_system_order[tso_index];
        if (fts.is_active(ts_index)) {
            sorted_active_ts_indices.push_back(ts_index);
            transition_system_label_ranks.push_back(vector<int>());
            vector<int> &label_ranks = transition_system_label_ranks.back();
            compute_label_ranks(fts, ts_index, label_ranks);
        }
    }

    pair<int, int> next_merge = get_next_dfp(fts, sorted_active_ts_indices);

    --remaining_merges;
    return next_merge;
}

pair<int, int> MergeDFP::get_next(FactoredTransitionSystem &fts,
                                  const vector<int> &ts_indices) {
    assert(initialized());
    assert(!done());

    /*
      Precompute a vector sorted_active_ts_indices which contains all exisiting
      transition systems in the given order and compute label ranks.
    */
    assert(!transition_system_order.empty());
    vector<int> sorted_active_ts_indices;
    for (size_t tso_index = 0; tso_index < transition_system_order.size(); ++tso_index) {
        int ts_index = transition_system_order[tso_index];
        for (int given_index : ts_indices) {
            if (ts_index == given_index) {
                assert(fts.is_active(ts_index));
                sorted_active_ts_indices.push_back(ts_index);
            }
        }
    }

    pair<int, int> next_merge = get_next_dfp(fts, sorted_active_ts_indices);

    --remaining_merges;
    return next_merge;
}

string MergeDFP::name() const {
    return "dfp";
}

void MergeDFP::compute_ts_order(const TaskProxy &task_proxy,
                                AtomicTSOrder atomic_ts_order,
                                ProductTSOrder product_ts_order,
                                bool atomic_before_product,
                                bool randomized_order,
                                vector<int> &ts_order) {
    int num_variables = task_proxy.get_variables().size();
    int max_transition_system_count = num_variables * 2 - 1;

    ts_order.reserve(max_transition_system_count);
    if (randomized_order) {
        for (int i = 0; i < max_transition_system_count; ++i) {
            ts_order.push_back(i);
        }
        g_rng.shuffle(ts_order);
    } else {

        // Compute the order in which atomic transition systems are considered
        vector<int> atomic_tso;
        for (int i = 0; i < num_variables; ++i) {
            atomic_tso.push_back(i);
        }
        if (atomic_ts_order == INVERSE) {
            reverse(atomic_tso.begin(), atomic_tso.end());
        } else if (atomic_ts_order == RANDOM1) {
            g_rng.shuffle(atomic_tso);
        }

        // Compute the order in which product transition systems are considered
        vector<int> product_tso;
        for (int i = num_variables; i < max_transition_system_count; ++i) {
            product_tso.push_back(i);
        }
        if (product_ts_order == NEW_TO_OLD) {
            reverse(product_tso.begin(), product_tso.end());
        } else if (product_ts_order == RANDOM2) {
            g_rng.shuffle(product_tso);
        }

        // Put the orders in the correct order
        if (atomic_before_product) {
            ts_order.insert(ts_order.end(),
                                           atomic_tso.begin(),
                                           atomic_tso.end());
            ts_order.insert(ts_order.end(),
                                           product_tso.begin(),
                                           product_tso.end());
        } else {
            ts_order.insert(ts_order.end(),
                                           product_tso.begin(),
                                           product_tso.end());
            ts_order.insert(ts_order.end(),
                                           atomic_tso.begin(),
                                           atomic_tso.end());
        }
    }
}

void MergeDFP::add_options_to_parser(OptionParser &parser, bool dfp_defaults) {
    vector<string> atomic_ts_order;
    atomic_ts_order.push_back("REGULAR");
    atomic_ts_order.push_back("INVERSE");
    atomic_ts_order.push_back("RANDOM");
    parser.add_enum_option("atomic_ts_order",
                           atomic_ts_order,
                           "order of atomic transition systems",
                           "REGULAR");
    vector<string> product_ts_order;
    product_ts_order.push_back("OLD_TO_NEW");
    product_ts_order.push_back("NEW_TO_OLD");
    product_ts_order.push_back("RANDOM");
    if (dfp_defaults) {
        parser.add_enum_option("product_ts_order",
                               product_ts_order,
                               "order of product transition systems",
                               "NEW_TO_OLD");
    } else {
        parser.add_enum_option("product_ts_order",
                               product_ts_order,
                               "order of product transition systems",
                               "OLD_TO_NEW");
    }
    if (dfp_defaults) {
        parser.add_option<bool>("atomic_before_product",
                                "atomic ts before product ts",
                                "false");
    } else {
        parser.add_option<bool>("atomic_before_product",
                                "atomic ts before product ts",
                                "true");
    }
    parser.add_option<bool>("randomized_order",
                            "globally randomized order",
                            "false");
}

static shared_ptr<MergeStrategy>_parse(OptionParser &parser) {
    MergeDFP::add_options_to_parser(parser);
    Options options = parser.parse();
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
    if (parser.dry_run())
        return nullptr;
    else
        return make_shared<MergeDFP>(options);
}

static PluginShared<MergeStrategy> _plugin("merge_dfp", _parse);
}
