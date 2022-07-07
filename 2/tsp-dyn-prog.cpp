#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <unordered_map>
#include <vector>

#include <boost/dynamic_bitset.hpp>
#include <boost/math/special_functions/binomial.hpp>

#define RESET "\033[0m"
#define RED "\033[31m"
#define GREEN "\033[32m"

// node/vertex in a graph
using node = uint32_t;

// a graph is represented as dense symmetric matrix
using Graph = std::vector<std::vector<double>>;

// infinite distance
constexpr double inf = std::numeric_limits<double>::max();

// Store values of gamma(S, u) for a fixed subset S a several nodes.
struct Gamma {
    Gamma() = default;

    // Constructs the Gamma struct for a subset S = {u} with the value
    // gamma(S, u) = gamma_val.
    Gamma(node u, double gamma_u) : nodes({u}), values({gamma_u}) {}

    // Nodes u in gamma(S, u)
    std::vector<node> nodes;

    // Values of gamma(S, u) of the nodes in `nodes`
    std::vector<double> values;

    // Inserts the value gamma(S, u) for a new node u
    void insert_gamma_value(node u, double gamma_val) {
        nodes.push_back(u);
        values.push_back(gamma_val);
    }

    // Iterates over node u and every value gamma(S, u)
    template <class Lambda>
    void iter_items(Lambda lambda) {
        for (uint32_t i = 0; i < nodes.size(); ++i)
            lambda(nodes[i], values[i]);
    }
};

// Hash map structure that stores a subset S and the values gamma(S, u)
using SubsetMap = std::unordered_map<boost::dynamic_bitset<>, Gamma>;

class TSP {
    // Input graph
    const Graph &G;

    // Number of nodes in the graph G
    const node n;

    // At index i stores every subset S with size `i`
    std::vector<SubsetMap> subsets;

    // Used to create every subset S
    std::vector<bool> permutation;
    boost::dynamic_bitset<> bitset;

    // Cheapest Hamiltonian path in G
    double path_cost;
    std::vector<node> path;

    // Creates every subset S with size `set_size` that also contain node 0.
    // Yields S as a vector of nodes (first parameter) and as a bitset (second parameter).
    template <class Lambda>
    void iter_subsets(uint32_t set_size, Lambda lambda) {
        assert(set_size > 1 && "Use this iterator for set_size > 1");
        std::fill(permutation.begin(), permutation.begin() + set_size - 1, true);
        std::fill(permutation.begin() + set_size - 1, permutation.end(), false);

        std::vector<node> subset;
        do {
            subset.clear();
            for (node u = 0; u < permutation.size(); ++u) {
                if (permutation[u])
                    subset.push_back(u + 1);
                bitset[u + 1] = permutation[u];
            }
            lambda(subset, bitset);
        } while (std::prev_permutation(permutation.begin(), permutation.end()));
    }

    // Stores the first subset {0} and gamma({0}, 0)
    void base_case() {
        // TODO: write your implementation here
        // A way to create a bitset with size n with all bits set to zero is
        // boost::dynamic_bitset<> my_bitset(n);

	    boost::dynamic_bitset<> myBitset(n);
	    myBitset[0] = 1;
	    Gamma g(0,0);
	    SubsetMap sm = { {myBitset, g} };
    	    subsets[1] = sm;
    }

    // Executes the i-th step of the dynamic programming (assuming that the previous steps are
    // completed).
    void dyn_prog(uint32_t i) {
        // TODO: write your implementation here

        // Avoid to access items of an unordered_map using the [] operator. To
        // access items that are expected to be in the map use 'at' because, if
        // the item is not found, it throws an exception. In this way it is
        // easier to find implementation errors.

	SubsetMap prev = subsets[i-1];
	SubsetMap sm;


	iter_subsets(i, [&](const std::vector<node> &S_vec, boost::dynamic_bitset<> S_bitset) {


			   Gamma g(0, inf);

			   for(node u: S_vec)
			   {
			   	if(u == 0) continue;



				S_bitset[u] = 0; //to find a Gamma Set where our destination u was not part of...
				Gamma prevGamma = prev.at(S_bitset);
				double findMin = inf;

				prevGamma.iter_items([&](node v, double v_gammaval){

					       if(v != u)
					       {
					       		double dist = v_gammaval + G[v][u];
					       		if ( dist < findMin) findMin = dist;
					       }
				});

				S_bitset[u] = 1; //Reset Bitset for further operations
				g.insert_gamma_value(u, findMin);

			   }

			   sm.insert( {S_bitset, g} );

	});
	subsets[i] = sm;

    }

    // Computes the Hamiltonian path with minimum cost and its cost.
    // Stores the results in the `path` and 'path_cost' member variable.
    void ham_path() {
        // TODO: write your implementation here

	path_cost = 0;
	path = std::vector<node>(n,0);

	node last = 0;


	boost::dynamic_bitset<> iterSet(n);
	iterSet.flip();

	for(int i = G.size(); i > 0 ; i--)
	{
		node minNode;
		double minVal = inf;

		Gamma g = subsets[i].at(iterSet);

		g.iter_items([&](node v, double v_gammaval){

				double dist = v_gammaval + G[v][last];
				if(dist < minVal)
				{
			       		minVal = dist;
					minNode = v;
				}
		});

		path[i-1] = minNode;

		if(last == 0) path_cost = minVal;

		iterSet[minNode] = 0;

		last = minNode;
	}

    }

    // No more tasks below this line.

    // Checks that the i-th dyn-prog step.
    void check_step(uint32_t i) const {
        uint32_t n_substes = boost::math::binomial_coefficient<double>(n - 1, i - 1);
        assert(subsets[i].size() == n_substes &&
               "The number of subsets with size i should be binomial(n - 1, i - 1)");
        for (const auto &S : subsets[i]) {
            assert(S.first.count() == i && "S must contain i elements including 0");
            assert(S.first[0] == 1 && "S must always contain node 0");
        }
    }

public:
    TSP(const Graph &G) : G(G), n(G.size()), bitset(n) {
        subsets.resize(n + 1);
        permutation.resize(n - 1);
        bitset[0] = 1;
        path.resize(n);
    }

    // Runs the algorithm
    void run() {
        base_case();

        for (uint32_t i = 2; i <= n; ++i) {
            dyn_prog(i);
            // Check correctness of step i
            check_step(i);
        }
        ham_path();
    }

    // Slow naive implementation (for debugging purposes)
    void naive() {
        std::vector<node> nodes;
        for (node u = 0; u < n; ++u)
            nodes.push_back(u);
        path_cost = inf;
        do {
            double cur_length = 0;
            for (uint32_t i = 0; i < n - 1; ++i)
                cur_length += G[nodes[i]][nodes[i + 1]];
            cur_length += G[nodes.front()][nodes.back()];
            if (cur_length < path_cost) {
                path_cost = cur_length;
                path = nodes;
            }
        } while (std::next_permutation(nodes.begin(), nodes.end()));
    }

    // Returns the length of the cheapest Hamiltonian path in G.
    double get_path_cost() const noexcept { return path_cost; }

    // Returns the cheapest Hamiltonian path in G.
    const std::vector<node> &get_path() const noexcept { return path; }
};

Graph read_graph(const std::string &path, uint32_t n) {
    std::ifstream file(path, std::ios::binary);
    Graph G(n);
    for (auto &row : G)
        row.resize(n);
    double val;
    node row = 0, col = 1;
    while (file.read(reinterpret_cast<char *>(&val), sizeof(double))) {
        G[row][col] = G[col][row] = static_cast<double>(val);
        if (col++ == n - 1) {
            ++row;
            col = row + 1;
        }
    }
    return G;
}

std::vector<node> read_sol(const std::string &path) {
    std::ifstream file(path, std::ios::binary);
    std::vector<node> sol;
    node next;
    while (file.read(reinterpret_cast<char *>(&next), sizeof(node)))
        sol.push_back(static_cast<node>(next));
    return sol;
}

template <class T>
void vector_to_file(const std::vector<T> &vec, const std::string &path) {
    std::ofstream outStream(path, std::ios::out | std::ios::binary);
    outStream.write(reinterpret_cast<const char *>(&vec[0]), vec.size() * sizeof(T));
    outStream.close();
}

std::string path_to_graph(uint32_t seed, uint32_t n) {
    return "./graphs/graph_" + std::to_string(n) + "_" + std::to_string(seed);
}

std::string path_to_sol(uint32_t seed, uint32_t n) { return path_to_graph(seed, n) + "_sol"; }

using EdgeList = std::vector<std::pair<node, node>>;

class Test {

    // Edge list of a cycle
    static EdgeList build_edge_list(const std::vector<node> &path) {
        const uint32_t n = path.size();
        EdgeList edge_list(n);
        for (uint32_t i = 0; i < n; ++i) {
            node u = path[i];
            node prev = (i == 0) ? path.back() : path[i - 1];
            node next = (i == n - 1) ? path.front() : path[i + 1];
            edge_list[u] = {std::min(prev, next), std::max(prev, next)};
        }
        return edge_list;
    }

    static double path_cost(const Graph &G, const std::vector<node> &path) {
        double cost = G[path.back()][path.front()];
        for (uint32_t i = 0; i < path.size() - 1; ++i)
            cost += G[path[i]][path[i + 1]];
        return cost;
    }

    static bool compare_edge_lists(const EdgeList &l1, const EdgeList &l2) {
        for (uint32_t i = 0; i < l1.size(); ++i)
            if (l1[i] != l2[i])
                return false;
        return true;
    }

public:
    static void test_graph(uint32_t n, uint32_t seed) {
        std::cout << "Test: " << n << " nodes, seed " << seed << ": ";
        const Graph G = read_graph(path_to_graph(seed, n), n);
        const auto sol_path = read_sol(path_to_sol(seed, n));
        const auto sol_edge_list = build_edge_list(sol_path);

        TSP tsp(G);
        tsp.run();
        const auto dyn_prog_path = tsp.get_path();
        const auto dyn_prog_edge_list = build_edge_list(dyn_prog_path);
        assert(dyn_prog_path.size() == n && "The path must visit all nodes in G once");
        if (!compare_edge_lists(sol_edge_list, dyn_prog_edge_list)) {
            std::cout << RED << "Failed\n" << RESET << "The computed path has cost ";
            std::cout << path_cost(G, dyn_prog_path) << " but a cheaper path has cost "
                      << path_cost(G, sol_path) << '\n';
        } else
            std::cout << GREEN << "Passed\n" << RESET;
    }
};

int main(int argc, char **argv) {
    for (uint32_t seed = 1; seed <= 5; ++seed)
        for (uint32_t n = 10; n <= 16; ++n)
            Test::test_graph(n, seed);

    return 0;
}
