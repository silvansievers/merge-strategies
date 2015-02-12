#include "graph_creator.h"

#include "../labels.h"
#include "../transition_system.h"

#include "../../bliss/defs.hh"

#include "../../globals.h"
#include "../../option_parser.h"
#include "../../timer.h"
#include "../../utilities.h"

#include <cassert>
#include <iostream>

using namespace std;

static void out_of_memory_handler() {
    throw bliss::BlissMemoryOut();
}

void add_automorphism(void* param, unsigned int, const unsigned int * automorphism) {
    GraphCreator *gc = (GraphCreator*) param;
    gc->create_symmetry_generator(automorphism);
}

GraphCreator::GraphCreator(const Options &options)
    : debug(options.get<bool>("debug_graph_creator")),
      stabilize_transition_systems(options.get<bool>("stabilize_transition_systems")),
      bliss_time_limit(options.get<double>("bliss_time_limit")),
      stop_after_no_symmetries(options.get<bool>("stop_after_no_symmetries")),
      num_identity_generators(0),
      bliss_limit_reached(false) {
}

GraphCreator::~GraphCreator() {
    delete_generators();
}

void GraphCreator::delete_generators() {
    for (size_t i = 0; i < symmetry_generators.size(); i++){
        if (symmetry_generators[i]) {
            delete symmetry_generators[i];
            symmetry_generators[i] = 0;
        }
    }

    symmetry_generators.clear();
    symmetry_generator_info.reset();
    num_identity_generators = 0;
}

void GraphCreator::create_symmetry_generator(const unsigned int *automorphism) {
    SymmetryGenerator* symmetry_generator = new SymmetryGenerator(symmetry_generator_info, automorphism);
    //Only if we have non-identity generator we need to save it into the list of generators
//  cout << "Checking for identity" << endl;
    if(!symmetry_generator->identity()) {
        symmetry_generators.push_back(symmetry_generator);
    } else {
        delete symmetry_generator;
        ++num_identity_generators;
//        if (num_identity_generators > stop_after_false_generated) {
//            cout << endl << "Problems with generating symmetry group! Too many false generators." << endl;
//            cout<<"Number of generators: 0"<<endl;
//            exit_with(EXIT_CRITICAL_ERROR);
//        }
    }
}

double GraphCreator::compute_generators(const vector<TransitionSystem *>& transition_systems) {
    // Find (non) transition system stabilized symmetries depending on the chosen option.

    cout << "Computing generators for " << (stabilize_transition_systems? "" : "non ")
         << "transition system stabilized symmetries" << endl;

    Timer timer;
    new_handler original_new_handler = set_new_handler(out_of_memory_handler);
    try {
        cout << "Creating the bliss graph..." << endl;
        bliss::Digraph bliss_graph = bliss::Digraph();
        bliss_graph.set_time_limit(bliss_time_limit);

        create_bliss_graph(transition_systems, bliss_graph);
    //    bliss_graph.set_splitting_heuristic(bliss::Digraph::shs_flm);
        bliss_graph.set_splitting_heuristic(bliss::Digraph::shs_fs);

        bliss::Stats stats1;
    //    FILE *f = fopen("bliss.log", "w");
    //    FILE *stats_file = fopen("bliss.stats", "w");
    //    bliss_graph.set_verbose_file(f);
    //    bliss_graph.set_verbose_level(10);

        cout << "Searching for automorphisms... " << timer << endl;
        bliss_graph.find_automorphisms(stats1,&(add_automorphism), this);
  //    stats1.print(stats_file);
  //    fclose(stats_file);
  //    fclose(f);

        cout << "Got " << symmetry_generators.size() << " group generators" << endl; //, time step: [t=" << g_timer << "]" << endl;
        cout << "Got " << num_identity_generators << " identity generators" << endl;
    } catch (bliss::BlissException &e) {
        e.dump();
        bliss_limit_reached = true;
    }
    set_new_handler(original_new_handler);

    if (stop_after_no_symmetries && symmetry_generators.empty()) {
        bliss_limit_reached = true;
    }

    cout << "Done computing symmetries: " << timer << endl;
    return timer();
}

void GraphCreator::create_bliss_graph(const vector<TransitionSystem *> &transition_systems,
                                      bliss::Digraph &bliss_graph) {
    int idx = 0;

    // Setting the number of transition systems
    symmetry_generator_info.num_transition_systems = transition_systems.size();;
    int node_color_added_val = 0;

    if (debug) {
        cout << "digraph pdg";
        cout << " {" << endl;
        cout << "    node [shape = none] start;" << endl;
    }

    // We start with one node for every transition system, later adding more nodes
    // for states and labels.
    int num_of_nodes = transition_systems.size();

    for (int abs_ind = 0; abs_ind < static_cast<int>(transition_systems.size()); ++abs_ind){
        // Add vertex for each transition system
        if (stabilize_transition_systems || transition_systems[abs_ind] == 0) {
            // Either the transition system is empty or all transition systems are stabilized.
            node_color_added_val++;
            // NOTE: we need to add an transition system vertex for every transition system, even the unused ones,
            // because we want to use the transition system indices as vertex-IDs and vertex IDs in a bliss graph
            // are numbered from 0 to n-1.
            // We further add an extra color for each empty transition system even when
            // when not stabilizing transition systems in order to ensure that no trivial
            // symmetries that map two empty transition systems to each other are found.
            idx = bliss_graph.add_vertex(TRANSITION_SYSTEM_VERTEX + node_color_added_val);
        } else {
          idx = bliss_graph.add_vertex(TRANSITION_SYSTEM_VERTEX);
        }
 //     cout << "Adding transition system vertex: " << idx << endl;
        assert(abs_ind == idx);
        if (debug) {
            cout << "    node" << idx << " [shape=circle, label=abs" << abs_ind << "]; // color: " << node_color_added_val << endl;
        }

        // Setting the indices for connections between abstract states and their transition systems
        symmetry_generator_info.dom_sum_by_var.push_back(num_of_nodes);

        int abs_states = 0;
        if (transition_systems[abs_ind])
            abs_states = transition_systems[abs_ind]->get_size();
        num_of_nodes += abs_states;
        for(int num_of_value = 0; num_of_value < abs_states; num_of_value++){
            symmetry_generator_info.var_by_val.push_back(abs_ind);
        }

    }
    // Setting the total number of abstract states and transition systems
    symmetry_generator_info.num_abs_and_states = num_of_nodes;

    // We need an arbitrary valid transition system to get access to the number of
    // labels and their costs (we do not have access to the labels object).
    const TransitionSystem *some_abs = 0;

    for (size_t abs_ind = 0; abs_ind < transition_systems.size(); abs_ind++){
        if (transition_systems[abs_ind] == 0)  //In case the transition system is empty
            continue;
        if (!some_abs)
            some_abs = transition_systems[abs_ind];

        int abs_states = transition_systems[abs_ind]->get_size();

        // Now add abs states for each transition system
        for(AbstractStateRef state = 0; state < abs_states; state++){

            idx = bliss_graph.add_vertex(ABSTRACT_STATE_VERTEX + node_color_added_val);
//          cout << "Added abstract state vertex: " << idx << " for transition system " << abs_ind << endl;

            // Edge from transition system nodes to all the transition system's states
            bliss_graph.add_edge(abs_ind, idx);

            if (debug) {
                cout << "    node" << idx << " [shape=circle, label=abs" << abs_ind << "_state" << state << "];" << endl;
                cout << "    node" << abs_ind << " -> node" << idx << ";" << endl;
            }
        }
    }

    // Now we add vertices for operators
    const Labels *labels = some_abs->get_labels();
    int num_labels = labels->get_size();
    for (int label_no = 0; label_no < num_labels; ++label_no){
        if (!labels->is_current_label(label_no))
            continue;
//      int label_op_by_cost = 3 * g_operators[op_no].get_cost();
        // Changed to one node per transition - two colors per operator
        int label_cost = 2 * labels->get_label_cost(label_no); // was label_op_by_cost

        // For each operator we have one label node
        int label_idx = bliss_graph.add_vertex(LABEL_VERTEX + label_cost + node_color_added_val);
//      cout << "Added label vertex: " << label_idx << " with color " << LABEL_VERTEX + label_op_by_cost <<" for operator " << op_no << endl;

        if (debug) {
            cout << "    node" << label_idx << " [shape=circle, label=label_no" << label_no /*<< "_" << g_operators[label_no].get_name()*/ << "];" << endl;
        }

        for (size_t abs_ind = 0; abs_ind < transition_systems.size(); abs_ind++){
            if (transition_systems[abs_ind] == 0)  //In case the transition system is empty
                continue;

            const std::vector<Transition>& transitions = transition_systems[abs_ind]->get_const_transitions_for_label(label_no);
            for (size_t i = 0; i < transitions.size(); ++i) {
                const Transition &trans = transitions[i];

                // For each transition we have a pair of nodes - pre and eff, both connected to their label node
                int pre_idx = bliss_graph.add_vertex(LABEL_VERTEX + label_cost + 1 + node_color_added_val);

//              cout << "Added pre vertex: " << pre_idx << " with color " << LABEL_VERTEX + label_op_by_cost + 1 <<" for operator " << op_no << " in transition system " << abs_ind << endl;
//                int eff_idx = g->add_vertex(LABEL_VERTEX + label_op_by_cost + 2);
                int eff_idx = pre_idx;
//              cout << "Added eff vertex: " << eff_idx << " with color " << LABEL_VERTEX + label_op_by_cost + 2 <<" for operator " << op_no << " in transition system " << abs_ind << endl;

                int src_idx = symmetry_generator_info.get_index_by_var_val_pair(abs_ind, trans.src);
                int target_idx = symmetry_generator_info.get_index_by_var_val_pair(abs_ind, trans.target);
                // Edges from abstract state source over pre=eff=transition-node to target
                bliss_graph.add_edge(src_idx, pre_idx);
                bliss_graph.add_edge(eff_idx, target_idx);

//                g->add_edge(pre_idx, eff_idx);

                // Edge from operator-label to transitions (pre=eff=transition-node) induced by that operator
                bliss_graph.add_edge(label_idx, pre_idx);

//                g->add_edge(label_idx, eff_idx);

                if (debug) {
                    cout << "    node" << pre_idx << " [shape=circle, label=transition];" << endl;
                    cout << "    node" << src_idx << " -> node" << pre_idx << ";" << endl;
                    cout << "    node" << eff_idx << " -> node" << target_idx << ";" << endl;
                    cout << "    node" << label_idx << " -> node" << pre_idx << ";" << endl;
                }
            }
        }
    }

    // Finally, adding the goal node
    idx = bliss_graph.add_vertex(GOAL_VERTEX + node_color_added_val);

    if (debug) {
        cout << "    node [shape = doublecircle] node" << idx << " [label = goal];" << endl;
    }

//  cout << "Added goal vertex: " << idx << " with color " << GOAL_VERTEX << endl;

    for (size_t abs_ind = 0; abs_ind < transition_systems.size(); abs_ind++){
        if (transition_systems[abs_ind] == 0)  //In case the transition system is empty
            continue;

        int abs_states = transition_systems[abs_ind]->get_size();

        for(AbstractStateRef state = 0; state < abs_states; state++){
            if (!transition_systems[abs_ind]->is_goal_state(state))
                continue;

            int val_idx = symmetry_generator_info.get_index_by_var_val_pair(abs_ind, state);

            // Edges from goal states to the goal node
            bliss_graph.add_edge(val_idx, idx);

            if (debug) {
                cout << "    node" << val_idx << " -> node" << idx << ";" << endl;
            }
        }
    }

    if (debug) {
        cout << "}" << endl;
    }
}
