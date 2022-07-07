#include <algorithm>
#include <cassert>
#include <cmath>
#include <random>
#include <stdexcept>
#include <unordered_set>
#include <vector>

#include "utils.hpp"

class TSP2OPT {

    const Graph &G;
    const node n;

    // To make the execution reproducible.
    std::mt19937_64 gen;

    // Current solution (i.e., a Hamiltonian cycle)
    std::vector<node> cycle;

    // Stores the indices 'i' and 'j' of two nodes in the cycle.
    // cost_decrement represents how much the cost of the cycle would decrease
    // if, instead of cycle[i] -> cycle[i + 1] and cycle[j] -> cycle[j + 1] we visit
    // cycle[i] -> cycle[j] and cycle[i + 1] -> cycle[j + 1].
    struct NodeSwap {
        uint32_t i = 0, j = 0;
        double cost_decrement = 0;

        NodeSwap() = default;

        NodeSwap(uint32_t i, uint32_t j, double cost_decrement)
            : i(i), j(j), cost_decrement(cost_decrement) {}

        void update(uint32_t new_i, uint32_t new_j, double new_cost_decrement) {
            i = new_i;
            j = new_j;
            cost_decrement = new_cost_decrement;
        }
    };

    // Index of the node after the one at index 'i' in the cycle.
    uint32_t next_index(uint32_t i) const noexcept { return i == n - 1 ? 0 : i + 1; }

    // Iterates over all pairs of nodes that can be exchanged.
    template <class Lambda>
    void iter_node_pairs(Lambda lambda) {
        for (uint32_t i = 0; i < n - 2; ++i)
            for (uint32_t j = i + 2; j < n - 1; ++j)
                lambda(cycle[i], i, cycle[j], j);
    }

    // Initializes the current solution uniformly at random.
    void random_init() {
        cycle.reserve(n);
        for(node i=0; i<n; i++) {
            cycle.push_back(i);
        }
        std::shuffle(cycle.begin(), cycle.end(), gen);
    }

    // Finds the swap of nodes that maximizes the decrement of the cost of the
    // Hamiltonian cycle and returns it.
    NodeSwap compute_best_swap() {
        NodeSwap result(0, 0, 0);
        iter_node_pairs([&](node v_i, uint32_t i, node v_j, uint32_t j) {
            double actual_distance = G[v_i][cycle.at(next_index(i))] + G[v_j][cycle.at(next_index(j))];
            double alt_distance = G[v_i][v_j] + G[cycle.at(next_index(i))][cycle.at(next_index(j))];
            double diff = actual_distance - alt_distance;
            if(diff > result.cost_decrement)
                result.update(i, j, diff);
        });
        return result;
    }

    // Swaps all the nodes from swap.i (excluded) to swap.j (included), i.e.,
    // reverses the cycle from index swap.i + 1 to swap.j.
    void apply_swap(const NodeSwap &swap) {
        std::reverse(std::begin(cycle) + next_index(swap.i),
                     std::begin(cycle) + next_index(swap.j));
    }

    // No more tasks below this line

    // Check that 'cycle' is a cycle.
    void check_cycle() const {
        assert(cycle.size() == n && "The cycle must include n nodes.");
        assert(std::unordered_set<node>(cycle.begin(), cycle.end()).size() == cycle.size() &&
               "The cycle must include all nodes only once.");
    }

public:
    TSP2OPT(const Graph &G, int seed = 42) : G(G), n(G.size()) {
        if (n < 4)
            throw std::runtime_error("The graph must have at least 4 vertices.");
        gen.seed(seed);
    }

    // Runs the 2-opt algorithm.
    void run() {
        random_init();
        check_cycle();
        do {
            const auto best_swap = compute_best_swap();
            if (best_swap.cost_decrement <= 0)
                return;
            apply_swap(best_swap);
            check_cycle();
        } while (true);
    }

    // Returns the current solution.
    const std::vector<node> &get_path() const noexcept { return cycle; }
};
