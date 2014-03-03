#include "scc.h"
#include <algorithm>
#include <vector>
using namespace std;

vector<vector<unsigned int> > SCC::get_result() {
    unsigned int node_count = graph.size();
    dfs_numbers.resize(node_count, -1);
    dfs_minima.resize(node_count, -1);
    stack_indices.resize(node_count, -1);
    stack.reserve(node_count);
    current_dfs_number = 0;

    for (unsigned int i = 0; i < node_count; i++)
        if (dfs_numbers[i] == -1)
            dfs(i);

    reverse(sccs.begin(), sccs.end());
    return sccs;
}

void SCC::dfs(unsigned int vertex) {
    unsigned int vertex_dfs_number = current_dfs_number++;
    dfs_numbers[vertex] = dfs_minima[vertex] = vertex_dfs_number;
    stack_indices[vertex] = stack.size();
    stack.push_back(vertex);

    const vector<unsigned int> &successors = graph[vertex];
    for (unsigned int i = 0; i < successors.size(); i++) {
        unsigned int succ = successors[i];
        unsigned int succ_dfs_number = dfs_numbers[succ];
        if (succ_dfs_number == -1) {
            dfs(succ);
            dfs_minima[vertex] = min(dfs_minima[vertex], dfs_minima[succ]);
        } else if (succ_dfs_number < vertex_dfs_number && stack_indices[succ] != -1) {
            dfs_minima[vertex] = min(dfs_minima[vertex], succ_dfs_number);
        }
    }

    if (dfs_minima[vertex] == vertex_dfs_number) {
        unsigned int stack_index = stack_indices[vertex];
        vector<unsigned int> scc;
        for (unsigned int i = stack_index; i < stack.size(); i++) {
            scc.push_back(stack[i]);
            stack_indices[stack[i]] = -1;
        }
        stack.erase(stack.begin() + stack_index, stack.end());
        sccs.push_back(scc);
    }
}
