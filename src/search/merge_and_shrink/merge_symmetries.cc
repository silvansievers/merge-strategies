#include "merge_symmetries.h"

#include "factored_transition_system.h"
#include "merge_selector.h"
#include "merge_tree.h"
#include "transition_system.h"

#include "symmetries/symmetry_generator.h"
#include "symmetries/symmetry_group.h"

#include "../options/options.h"

#include <algorithm>
#include <iomanip>

using namespace std;

namespace merge_and_shrink {
MergeSymmetries::MergeSymmetries(
    const FactoredTransitionSystem &fts,
    const options::Options &options,
    int num_merges,
    unique_ptr<MergeTree> merge_tree,
    shared_ptr<MergeSelector> merge_selector,
    bool tree_is_miasm)
    : MergeStrategy(fts),
      num_merges(num_merges),
      symmetries_for_merging(SymmetriesForMerging(options.get_enum("symmetries_for_merging"))),
      internal_merging(InternalMerging(options.get_enum("internal_merging"))),
      max_bliss_iterations(options.get<int>("max_bliss_iterations")),
      bliss_call_time_limit(options.get<int>("bliss_call_time_limit")),
      tree_is_miasm(tree_is_miasm),
      bliss_remaining_time_budget(options.get<int>("bliss_total_time_budget")),
      merge_tree(move(merge_tree)),
      merge_selector(merge_selector),
      iteration_counter(0),
      bliss_limit_reached(false),
      pure_fallback_strategy(true) {
    symmetry_group = utils::make_unique_ptr<SymmetryGroup>(
        options.get<bool>("debug_graph_creator"),
        options.get<bool>("stabilize_transition_systems"));
}

void MergeSymmetries::dump_statistics() {
    if (!bliss_times.empty()) {
        sort(bliss_times.begin(), bliss_times.end());
        double summed_up_bliss_times = 0;
        size_t total_bliss_calls = bliss_times.size();
        cout << "Total bliss calls: " << total_bliss_calls << endl;
        for (size_t i = 0; i < total_bliss_calls; ++i) {
            summed_up_bliss_times += bliss_times[i];
        }
        double average_bliss_time = summed_up_bliss_times
                / (double) total_bliss_calls;
        cout << setprecision(5) << fixed << "Average bliss time: " << average_bliss_time << endl;
        double median_bliss_time;
        if (total_bliss_calls % 2 == 0) {
            size_t lower_median_index = (total_bliss_calls - 1) / 2;
            size_t higher_median_index = lower_median_index + 1;
            median_bliss_time = (bliss_times[lower_median_index]
                                 + bliss_times[higher_median_index])
                    / 2.0;
        } else {
            size_t median_index = total_bliss_calls / 2;
            median_bliss_time = bliss_times[median_index];
        }
        cout << setprecision(5) << fixed << "Median bliss time: " << median_bliss_time << endl;
    }
}

void MergeSymmetries::determine_merge_order() {
    int chosen_generator_for_merging = -1;
    int smallest_generator_affected_transition_systems_size = numeric_limits<int>::max();
    int largest_generator_affected_transition_systems_size = 0;

    /*
      Go over the generators and classify them into atomic, local or general
      ones. Also compute the generator that affects or maps the smallest or
      largest number of transition systems (depending on the chosen options).
    */
    for (int generator_index = 0; generator_index < symmetry_group->get_num_generators(); ++generator_index) {
        const SymmetryGenerator *generator = symmetry_group->get_generator(generator_index);
        const vector<int> &overall_affected_transition_systems =
            generator->get_overall_affected_transition_systems();

        int number_overall_affected_transition_systems = overall_affected_transition_systems.size();
        if (number_overall_affected_transition_systems < 1) {
            cerr << "Something is wrong! The generator is the identity generator." << endl;
            utils::exit_with(utils::ExitCode::CRITICAL_ERROR);
        }
        if (number_overall_affected_transition_systems > 1) {
            if (symmetries_for_merging == SMALLEST
                    && number_overall_affected_transition_systems
                    < smallest_generator_affected_transition_systems_size) {
                smallest_generator_affected_transition_systems_size
                        = number_overall_affected_transition_systems;
                chosen_generator_for_merging = generator_index;
            } else if (symmetries_for_merging == LARGEST
                          && number_overall_affected_transition_systems
                          > largest_generator_affected_transition_systems_size) {
                largest_generator_affected_transition_systems_size
                        = number_overall_affected_transition_systems;
                chosen_generator_for_merging = generator_index;
            }
        }

        // Dump generator properties -- warning! lots of output!
//        const vector<int> &internally_affected_transition_systems =
//            generator->get_internally_affected_transition_systems();
//        cout << "Generator " << generator_index << endl;
//        for (size_t i = 0; i < mapped_transition_systems.size(); ++i) {
//            int abs_index = mapped_transition_systems[i];
//            int to_index = generator->get_value(abs_index);
//            cout << transition_systems[abs_index]->tag() << " mapped to " <<
//                    transition_systems[to_index]->tag();
//            if (generator->internally_affects(abs_index))
//                cout << " (and also internally affected)";
//            cout << endl;
//        }
//        for (size_t i = 0; i < internally_affected_transition_systems.size(); ++i) {
//            int abs_index = internally_affected_transition_systems[i];
//            assert(!generator->maps(abs_index));
//            cout << transition_systems[abs_index]->tag() << " internally affected" << endl;
//        }
    }

    /*
      Determine a merge strategy as follows:
      If we merge towards local symmetries, we only need to merge all mapped
      transition systems. They are always merged according to the cycle(s)
      of those mapped transition systems. If we merge towards atomic symmetries,
      it depends on whether we want to merge all the affected transition systems
      linearly or non-linearly: if linearly, we put them in the order given by
      fast downward. If non-linearly, we first merge all cycles (as when merging
      for local symmetries), and then linearly the products of these merges
      and the remaining other (internally) affected transition systems, again
      as given by fast downward.
    */
    if (symmetries_for_merging != NO_MERGING && chosen_generator_for_merging != -1) {
        vector<vector<int> > cycles;
        vector<int> merge_linear_transition_systems;
        const SymmetryGenerator *generator = symmetry_group->get_generator(chosen_generator_for_merging);

        // Always include all mapped transition systems
        if (internal_merging == NON_LINEAR) {
            /*
              If the internal merge strategy is non linear or we only want
              to merge every cycle (non linearly), we need to
              compute the actual cycles of transition system mappings.
            */
            generator->compute_cycles(cycles);
        } else if (internal_merging == LINEAR) {
            /*
              If the internal merge strategy is linear, we simply collect
              all mapped transition systems (i.e. the same transition systems
              as above, but we do not compute cycle information)
            */
            const vector<int> &mapped_transition_systems =
                generator->get_mapped_transition_systems();
            merge_linear_transition_systems.insert(merge_linear_transition_systems.end(),
                                                   mapped_transition_systems.begin(),
                                                   mapped_transition_systems.end());
        }

        /*
          We need to include the internally affected transition systems (to be
          merged linearly after the mapped ones).
        */
        const vector<int> &internally_affected_transition_systems =
             generator->get_internally_affected_transition_systems();
        merge_linear_transition_systems.insert(merge_linear_transition_systems.end(),
                                               internally_affected_transition_systems.begin(),
                                               internally_affected_transition_systems.end());

        /*
          Here we compute the actual merge tree: if cycles is non-empty, we
          start by merging all mapped transition systems cycle-wise. At the
          same time, we need to remember the indices of the transition systems
          resulting from those merges in case we want to merge those afterwards
          (when merging for atomic symmetries).
        */
        assert(merge_order.empty());
        int number_of_transition_systems = fts.get_size();
        int number_of_merges = 0;
        vector<int> merge_linear_indices;
        for (size_t cycle_no = 0; cycle_no < cycles.size(); ++cycle_no) {
            const vector<int> &cycle = cycles[cycle_no];
            size_t abs_index_1 = cycle[0];
            for (size_t i = 1; i < cycle.size(); ++i) {
                size_t abs_index_2 = cycle[i];
                merge_order.push_back(make_pair(abs_index_1, abs_index_2));
                abs_index_1 = number_of_transition_systems + number_of_merges;
                ++number_of_merges;
            }
            /*
              number_of_transition_systems + number_of_merges always is the *next*
              position where a new merged transition system will be stored at.
              here, we need the *last* position where the transition system
              resulting from merging the cycle was stored, hence the -1.
            */
            merge_linear_indices.push_back(number_of_transition_systems + number_of_merges - 1);
        }

        /*
          Here we include the transition systems that need to be merged linearly
          in the merge tree. Those are the internally affected ones and the
          products of the merged cycles if applicable
          (merge_linear_transition_systems).
        */
        merge_linear_indices.insert(merge_linear_indices.end(),
                                    merge_linear_transition_systems.begin(),
                                    merge_linear_transition_systems.end());

        size_t abs_index_1 = merge_linear_indices[0];
        for (size_t i = 1; i < merge_linear_indices.size(); ++i) {
            size_t abs_index_2 = merge_linear_indices[i];
            merge_order.push_back(make_pair(abs_index_1, abs_index_2));
            abs_index_1 = number_of_transition_systems + number_of_merges;
            ++number_of_merges;
        }

        cout << "Chosen internal merge order: " << endl;
        for (size_t i = 0; i < merge_order.size(); ++i) {
            cout << merge_order[i].first << ", " << merge_order[i].second << endl;
        }
    }
}

pair<int, int> MergeSymmetries::get_next() {
    ++iteration_counter;

    if (!bliss_limit_reached && iteration_counter <= max_bliss_iterations && merge_order.empty()) {
        double time_limit = 0.0;
        if (bliss_remaining_time_budget) {
            // we round everything down to the full second
            time_limit = bliss_remaining_time_budget;
        } else if (bliss_call_time_limit) {
            time_limit = bliss_call_time_limit;
        }
        cout << "Setting bliss time limit to " << time_limit << endl;
        double bliss_time = symmetry_group->find_symmetries(fts, time_limit);
        if (symmetry_group->is_bliss_limit_reached()) {
            bliss_limit_reached = true;
            // We do not use the possibly incomplete generators (TODO: this
            // probably does not make sense, but we didn't use them before
            // so we stick to this policy).
        } else {
            // Compute merge order.
            determine_merge_order();
        }
        symmetry_group->reset();
        if (pure_fallback_strategy && (!merge_order.empty())) {
            pure_fallback_strategy = false;
            cout << "not pure fallback strategy anymore" << endl;
        }
        bliss_times.push_back(bliss_time);
        if (bliss_remaining_time_budget) {
            bliss_remaining_time_budget -= bliss_time;
            if (bliss_remaining_time_budget <= 0) {
                assert(bliss_limit_reached);
                // in case of measurement inaccuracies, set bliss limit
                // reached to true in any case.
                bliss_limit_reached = true;
            }
            cout << "Remaining bliss time budget " << bliss_remaining_time_budget << endl;
        }
    }

    if (iteration_counter == num_merges) { // TODO: check correctness
        dump_statistics();
    }

    int next_merge_index = fts.get_size();

    // No symmetries merge order -- use fallback strategy.
    if (merge_order.empty()) {
        if (merge_selector) {
            return merge_selector->select_merge(fts);
        } else {
            assert(merge_tree);
            return merge_tree->get_next_merge(next_merge_index);
        }
    }

    pair<int, int> next_merge = merge_order.front();
    if (!fts.is_active(next_merge.first) || !fts.is_active(next_merge.second)) {
        cerr << "Problem with the merge strategy: invalid indices" << endl;
        utils::exit_with(utils::ExitCode::CRITICAL_ERROR);
    }
    merge_order.erase(merge_order.begin());

    // Update the merge tree, if present.
    if (merge_tree) {
        merge_tree->update(
            next_merge,
            next_merge_index,
            tree_is_miasm);
    }

    return next_merge;
}
}
