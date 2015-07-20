#ifndef MIASM_MERGE_AND_SHRINK_ABSTRACTION_H
#define MIASM_MERGE_AND_SHRINK_ABSTRACTION_H

#include "types.h"

#include "../task_proxy.h"

#include <memory>
#include <string>
#include <vector>
#include <map>

class Labels;
class MergeStrategy;
class Options;
class ShrinkStrategy;
class TransitionSystem;
class VarSetInfoRegistry;

class MiasmAbstraction {
    const std::shared_ptr<AbstractTask> task;
    TaskProxy task_proxy;
    MergeStrategy *const merge_strategy;
    ShrinkStrategy *const shrink_strategy;
    /** @brief the transition labels */
    Labels *labels;
public:
    MiasmAbstraction(const Options &opts);
    virtual ~MiasmAbstraction();
    static std::string option_key();
    static std::string plugin_key();



    std::map<mst::var_set_t, TransitionSystem *> cache;
    void release_cache();
    void release_cache(const mst::var_set_t &var_set);

    TransitionSystem *build_transition_system(
        const mst::var_set_t &G,
        std::vector<mst::var_set_t> &newly_built,
        const VarSetInfoRegistry &vsir);

};

#endif // MIASM_MERGE_AND_SHRINK_ABSTRACTION_H
