#ifndef TASKS_ROOT_TASK_H
#define TASKS_ROOT_TASK_H

#include "../abstract_task.h"

#include <set>
#include <vector>

namespace tasks {
extern std::shared_ptr<AbstractTask> g_root_task;
extern void read_root_task(std::istream &in);
extern std::vector<std::vector<bool>> g_mutex_var_pairs;
extern std::vector<std::vector<std::set<FactPair>>> g_mutexes;
}
#endif
