#include "merge_tree_factory_manual.h"

#include "merge_tree.h"

#include "../task_proxy.h"

#include "../plugins/plugin.h"
#include "../plugins/options.h"

#include "../utils/logging.h"
#include "../utils/system.h"

#include <cassert>
#include <iostream>
#include <map>
#include <sstream>

using namespace std;
using utils::ExitCode;

namespace merge_and_shrink {
MergeTreeFactoryManual::MergeTreeFactoryManual(
    const plugins::Options &options)
    : MergeTreeFactory(options) {
    if (options.contains("merge_order_list")) {
        merge_order_list = options.get_list<vector<int>>("merge_order_list");
    }
    if (options.contains("merge_order_tree_string")) {
        merge_order_tree_string = options.get<string>("merge_order_tree_string");
    }
}

unique_ptr<MergeTree> MergeTreeFactoryManual::compute_merge_tree(
    const TaskProxy &task_proxy) {
    int num_vars = task_proxy.get_variables().size();
    MergeTreeNode *root = nullptr;
    if (!merge_order_list.empty()) {
        int num_merges = num_vars - 1;
        if (static_cast<int>(merge_order_list.size()) != num_merges) {
            cerr << "Number of merges in the given task: "
                 << num_merges << endl;
            cerr << "Number of merges in the specified merge order: "
                 << merge_order_list.size() << endl;
            cerr << "Invalid size of merge order" << endl;
            utils::exit_with(ExitCode::SEARCH_INPUT_ERROR);
        }
        map<int, MergeTreeNode *> index_to_tree;
        for (int atomic_ts_index = 0; atomic_ts_index < num_vars; ++atomic_ts_index) {
            index_to_tree[atomic_ts_index] = new MergeTreeNode(atomic_ts_index);
        }
        int next_ts_index = num_vars;
        for (const vector<int> &merge : merge_order_list) {
            assert(merge.size() == 2);
            int ts_index1 = merge[0];
            int ts_index2 = merge[1];
            index_to_tree[next_ts_index] =
                new MergeTreeNode(index_to_tree[ts_index1], index_to_tree[ts_index2]);
            ++next_ts_index;
        }
        root = index_to_tree[next_ts_index - 1];
    } else {
        assert(!merge_order_tree_string.empty());
        // clamp the first opening and the last closing bracket away
        root = new MergeTreeNode(merge_order_tree_string.substr(
                                     1, merge_order_tree_string.size() - 2));
    }

    return utils::make_unique_ptr<MergeTree>(root, rng, update_option);
}

string MergeTreeFactoryManual::name() const {
    return "manual";
}

void MergeTreeFactoryManual::dump_tree_specific_options(utils::LogProxy &log) const {
    if (!merge_order_list.empty()) {
        log << "given merge order, as list: " << merge_order_list << endl;
    } else {
        assert(!merge_order_tree_string.empty());
        log << "given merge order, as string: " << merge_order_tree_string << endl;
    }
}

/*
class MergeTreeFactoryManualFeature : public plugins::TypedFeature<MergeTreeFactory, MergeTreeFactoryManual> {
public:
    MergeTreeFactoryManualFeature() : TypedFeature("manual") {
        document_title("Manuel merge trees");
        document_synopsis(
            "Manually specify a merge tree either as a list of merges or a "
            "specific string describing a tree.");
        add_list_option<vector<int>>(
            "merge_order_list",
            "merge order as list. NOTE/TODO: the resulting merge tree cannot be"
            "guaranteed to be processed in the exact same order as the list.",
            plugins::ArgumentInfo::NO_DEFAULT);
        add_option<string>(
            "merge_order_tree_string",
            "merge tree, specified as a string of the following form: "
            "xx<X>yxx<Y>yx<Z>yyy, where x and y are opening and closing brackets, "
            "respectively, and X, Y and Z are integer values denoting variables "
            "(i.e. indices of atomic transition systems)",
            plugins::ArgumentInfo::NO_DEFAULT);

        MergeTreeFactory::add_options_to_feature(*this);
        MergeTreeFactoryManual::add_options_to_feature(*this);
    }

    virtual shared_ptr<MergeTreeFactoryManual> create_component(const plugins::Options &options, const utils::Context &context) const override {
        if (options.contains("merge_order_list") && options.contains("merge_order_tree_string")) {
            context.error("Specifying a merge order and a merge tree is not possible!");
        } else if (!options.contains("merge_order_list") && !options.contains("merge_order_tree_string")) {
            context.error("Neither a merge order nor a merge tree was specified!");
        }
        if (options.contains("merge_order_list")) {
            vector<vector<int>> merge_order = options.get_list<vector<int>>("merge_order_list");
            if (merge_order.empty()) {
                context.error("Got empty merge order, aborting");
            }
            for (const vector<int> &pair : merge_order) {
                if (pair.size() != 2) {
                    stringstream error_msg;
                    error_msg << "Every element in the list merge_order_list must "
                                 "contain exactly two elements!" << endl;
                    error_msg << pair << endl;
                    context.error(error_msg.str());
                }
            }
        }
        return make_shared<MergeTreeFactoryManual>(options);
    }
};

static plugins::FeaturePlugin<MergeTreeFactoryManualFeature> _plugin;
*/
}
