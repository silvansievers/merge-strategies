#include "merge_strategy_factory_predefined.h"

#include "factored_transition_system.h"
#include "merge_predefined.h"

#include "../options/option_parser.h"
#include "../options/options.h"
#include "../options/plugin.h"

#include "../utils/logging.h"
#include "../utils/system.h"

#include <cassert>
#include <iostream>

using namespace std;
using utils::ExitCode;

namespace merge_and_shrink {
class BinaryTree {
    int index;
    BinaryTree *left_child;
    BinaryTree *right_child;

    bool is_leaf() const {
        return !left_child && !right_child;
    }
    bool has_two_leaf_children() const {
        return left_child && right_child && left_child->is_leaf() && right_child->is_leaf();
    }
public:
    BinaryTree(string tree_string)
        : index(-1),
          left_child(nullptr),
          right_child(nullptr) {
        assert(!tree_string.empty());
//        cout << "building tree for string " << tree_string << endl;
        char &c = tree_string[0];
        if (c == 'y') {
            cerr << "ill-specified string" << endl;
            utils::exit_with(ExitCode::INPUT_ERROR);
        } else if (c == 'x') {
            assert(tree_string.size() > 5); // need to have at least two subtrees
            int parentheses_counter = 1;
            string::size_type index = 1;
            while (parentheses_counter != 0) {
                if (tree_string[index] == 'x') {
                    ++parentheses_counter;
                } else if (tree_string[index] == 'y') {
                    --parentheses_counter;
                }
                ++index;
            }
            string left_sub_tree_string = tree_string.substr(1, index - 2);
            left_child = new BinaryTree(left_sub_tree_string);

            assert(tree_string[index] == 'x');
            parentheses_counter = 1;
            ++index;
            int right_index = index;
            while (parentheses_counter != 0) {
                if (tree_string[index] == 'x') {
                    ++parentheses_counter;
                } else if (tree_string[index] == 'y') {
                    --parentheses_counter;
                }
                ++index;
            }
            assert(index == tree_string.size());
            string right_sub_tree_string = tree_string.substr(right_index,
                                                              index - right_index - 1);
            right_child = new BinaryTree(right_sub_tree_string);
        } else {
            string to_be_index = tree_string.substr(0, tree_string.size() - 2);
//            cout << "to be index: " << to_be_index << endl;
            index = stoi(to_be_index);
        }
    }

    ~BinaryTree() {
        delete left_child;
        delete right_child;
    }

    void get_next_merge(int new_index, int &next_index1, int &next_index2) {
        if (next_index1 != -1 && next_index2 != -1) {
            return;
        }
        if (has_two_leaf_children()) {
            assert(index == -1);
            index = new_index;
            next_index1 = left_child->index;
            next_index2 = right_child->index;
            delete left_child;
            delete right_child;
            left_child = nullptr;
            right_child = nullptr;
        } else {
            // TODO: this seems to be incorrect. If the left subtree "returns"
            // a next merge pair, we should be done.
            if (left_child) {
                left_child->get_next_merge(new_index, next_index1, next_index2);
            }
            if (right_child) {
                right_child->get_next_merge(new_index, next_index1, next_index2);
            }
        }
    }

    int compute_size() const {
        if (is_leaf()) {
            return 0;
        } else {
            int number_of_internal_nodes = 1; // count the node itself
            if (left_child) {
                number_of_internal_nodes += left_child->compute_size();
            }
            if (right_child) {
                number_of_internal_nodes += right_child->compute_size();
            }
            return number_of_internal_nodes;
        }
    }

    void postorder(int indentation) const {
        if (left_child) {
            left_child->postorder(indentation + 1);
        }
        if (right_child) {
            right_child->postorder(indentation + 1);
        }
        for (int i = 0; i < indentation; ++i) {
            cout << "  ";
        }
        cout << index << endl;
    }
};

MergeStrategyFactoryPredefined::MergeStrategyFactoryPredefined(
    const options::Options &options)
    : MergeStrategyFactory() {
    if (options.contains("merge_order")) {
        merge_order = options.get_list<vector<int>>("merge_order");
        root = nullptr;
    } else if (options.contains("merge_tree")) {
        string tree_string = options.get<string>("merge_tree");
        // clamp the first opening and the last closing bracket away
        root = new BinaryTree(tree_string.substr(1, tree_string.size() - 2));
    }
}

MergeStrategyFactoryPredefined::~MergeStrategyFactoryPredefined() {
    delete root;
}

unique_ptr<MergeStrategy> MergeStrategyFactoryPredefined::compute_merge_strategy(
    shared_ptr<AbstractTask>,
    FactoredTransitionSystem &fts) {
    int num_ts = fts.get_size();
    int num_merges = num_ts - 1;
    if (!merge_order.empty() &&
            static_cast<int>(merge_order.size()) != num_merges) {
        cout << "Number of merges in the given task: "
             << num_merges << endl;
        cout << "Number of merges in the specified merge order: "
             << merge_order.size() << endl;
        cerr << "Invalid size of merge order" << endl;
        utils::exit_with(ExitCode::INPUT_ERROR);
    }
    if (root && root->compute_size() != num_merges) {
        cout << "Number of merges in the given task: "
             << num_merges << endl;
        cout << "Number of merges in the specified merge tree: "
             << root->compute_size() << endl;
        cerr << "Invalid size of merge tree!" << endl;
        utils::exit_with(ExitCode::INPUT_ERROR);
    }

    if (root) {
        int tree_size = root->compute_size();
        int new_index = num_ts;
        merge_order.reserve(tree_size);
        while (tree_size > 0) {
            int next_index1 = -1;
            int next_index2 = -1;
            root->get_next_merge(new_index, next_index1, next_index2);
            vector<int> next_merge;
            next_merge.push_back(next_index1);
            next_merge.push_back(next_index2);
            merge_order.push_back(next_merge);
            tree_size = root->compute_size();
            ++new_index;
        }
    }

    return utils::make_unique_ptr<MergePredefined>(fts, move(merge_order));
}

void MergeStrategyFactoryPredefined::dump_strategy_specific_options() const {
    if (!merge_order.empty()) {
        cout << "merge order: " << merge_order << endl;
    } else if (root) {
        cout << "merge tree: " << endl;
        root->postorder(0);
    }
}

string MergeStrategyFactoryPredefined::name() const {
    return "predefined";
}

static shared_ptr<MergeStrategyFactory>_parse(options::OptionParser &parser) {
    parser.document_synopsis(
        "Random merge strategy.",
        "This merge strategy randomly selects the two next transition systems"
        "to merge.");
    parser.add_list_option<vector<int>>(
        "merge_order",
        "predefined merge order",
        options::OptionParser::NONE);
    parser.add_option<string>(
        "merge_tree",
        "merge tree, specified as a string of the following form: "
        "xx<X>yxx<Y>yx<Z>yyy, where x and y are opening and closing brackets, "
        "respectively, and X, Y and Z are integer values denoting variables "
        "(i.e. indices of atomic transition systems)",
        options::OptionParser::NONE);

    options::Options options = parser.parse();
    if (parser.dry_run()) {
        if (options.contains("merge_order") && options.contains("merge_tree")) {
            cerr << "Specifying a merge order and a merge tree is not possible!" << endl;
            utils::exit_with(ExitCode::INPUT_ERROR);
        } else if (!options.contains("merge_order") && !options.contains("merge_tree")) {
            cerr << "Neither a merge order nor a merge tree was specified!" << endl;
            utils::exit_with(ExitCode::INPUT_ERROR);
        }
        if (options.contains("merge_order")) {
            vector<vector<int>> merge_order = options.get_list<vector<int>>("merge_order");
            if (merge_order.empty()) {
                cerr << "Got empty merge order, aborting" << endl;
                utils::exit_with(ExitCode::INPUT_ERROR);
            }
            for (const vector<int> &pair : merge_order) {
                if (pair.size() != 2) {
                    cerr << "Every element in the list merge_order must contain "
                        "exactly two elements!" << endl;
                    cout << pair << endl;
                    utils::exit_with(ExitCode::INPUT_ERROR);
                }
            }
        }
        return nullptr;
    } else {
        return make_shared<MergeStrategyFactoryPredefined>(options);
    }
}

static options::PluginShared<MergeStrategyFactory> _plugin("merge_predefined", _parse);
}
