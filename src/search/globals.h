#ifndef GLOBALS_H
#define GLOBALS_H

#include <iosfwd>
#include <memory>
#include <set>
#include <string>
#include <vector>

class AbstractTask;
class Axiom;
class AxiomEvaluator;
class CausalGraph;
struct FactPair;
class GlobalOperator;
class GlobalState;
class StateRegistry;
class SuccessorGenerator;

namespace int_packer {
class IntPacker;
}

namespace utils {
struct Log;
class RandomNumberGenerator;
}

bool test_goal(const GlobalState &state);
/*
  Set generates_multiple_plan_files to true if the planner can find more than
  one plan and should number the plans as FILENAME.1, ..., FILENAME.n.
*/
void save_plan(const std::vector<const GlobalOperator *> &plan,
               bool generates_multiple_plan_files = false);
int calculate_plan_cost(const std::vector<const GlobalOperator *> &plan);

void read_everything(std::istream &in);
void dump_everything();

// The following six functions are deprecated. Use task_tools.h instead.
bool is_unit_cost();
bool has_axioms();
void verify_no_axioms();
bool has_conditional_effects();
void verify_no_conditional_effects();
void verify_no_axioms_no_conditional_effects();

void check_magic(std::istream &in, std::string magic);

bool are_mutex(const FactPair &a, const FactPair &b);

extern bool g_use_metric;
extern int g_min_action_cost;
extern int g_max_action_cost;

// TODO: The following five belong into a new Variable class.
extern std::vector<std::string> g_variable_name;
extern std::vector<int> g_variable_domain;
extern std::vector<std::vector<std::string>> g_fact_names;
extern std::vector<int> g_axiom_layers;
extern std::vector<int> g_default_axiom_values;

extern int_packer::IntPacker *g_state_packer;
// This vector holds the initial values *before* the axioms have been evaluated.
// Use a state registry to obtain the real initial state.
extern std::vector<int> g_initial_state_data;
extern std::vector<std::pair<int, int>> g_goal;

extern std::vector<GlobalOperator> g_operators;
extern std::vector<GlobalOperator> g_axioms;
extern AxiomEvaluator *g_axiom_evaluator;
extern SuccessorGenerator *g_successor_generator;
extern std::string g_plan_filename;
extern int g_num_previously_generated_plans;
extern bool g_is_part_of_anytime_portfolio;

extern const std::shared_ptr<AbstractTask> g_root_task();
extern std::vector<std::vector<bool>> g_mutex_var_pairs;
extern std::vector<std::vector<std::set<FactPair>>>
g_inconsistent_facts;

extern utils::Log g_log;

#endif
