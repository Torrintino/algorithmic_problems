#ifndef MAX_CUT_HPP_
#define MAX_CUT_HPP_

#include <cassert>
#include <iostream>
#include <random>

#include <networkit/auxiliary/BucketPQ.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/graph/GraphTools.hpp>

// If needed, add further includes

#define RESET "\033[0m"
#define RED "\033[31m"
#define GREEN "\033[32m"

using NetworKit::Graph;
using NetworKit::node;

class MaxCutGreedy final {
    const Graph *G;
    const int seed;
    std::vector<bool> partition;

    // Returns a random cut as a set of nodes picked at random.
    void generate_random_cut() {
       
	srand(seed);	

	G->forNodes([&](node u){

			if( ((float)rand()/RAND_MAX) > 0.5) partition[u] = true;		
	});
			
    }

    // Returns the key of node u in the min-priority queue.
    int64_t key_of_node(node u) const {
	

	 int64_t key = 0;
	 G->forNeighborsOf(u, [&](node neigh){
			
			 if(partition[neigh] == partition[u]) key--;
			 else key++;

	});

	return key;
    }

public:
    MaxCutGreedy(const Graph *G, int seed = 42)
        : G(G), seed(seed), partition(G->upperNodeIdBound()) {}

    std::vector<node> run_max_cut() {
        
	std::vector<node> A;

	generate_random_cut();
	
	Aux::BucketPQ prioQ(G->numberOfNodes(), -1 * NetworKit::GraphTools::maxDegree(*G), NetworKit::GraphTools::maxDegree(*G));

	G->forNodes([&](node u){
			prioQ.insert(key_of_node(u), u);
			});

	std::pair<int64_t,node> min = prioQ.extractMin();
	while(min.first < 0)
	{	
		partition[min.second] = !partition[min.second];

		G->forNeighborsOf(min.second, [&](node neigh){
				prioQ.changeKey(key_of_node(neigh), neigh);
				});

		min = prioQ.extractMin();
	}
	
	G->forNodes([&](node u){
		if(partition[u] == true) A.push_back(u);
		});
	
	return A;
    }

    // No further tasks below this line
};

class TestMaxCutGreedy final {
public:
    static void run_test(const NetworKit::Graph *G) {
        std::cout << "Testing graph with " << G->numberOfNodes() << " nodes... ";
        std::vector<bool> inCut(G->upperNodeIdBound());

        MaxCutGreedy mc(G);
        const auto partition = mc.run_max_cut();
        for (auto u : partition)
            inCut[u] = true;

        uint32_t cut_size = 0;
        int64_t in_same_partition, in_other_partition;
        for (node u : G->nodeRange()) {
            in_same_partition = in_other_partition = 0;
            G->forNeighborsOf(u, [&](node v) {
                if (inCut[u] == inCut[v])
                    ++in_same_partition;
                else
                    ++in_other_partition;
            });

            if (in_same_partition > in_other_partition) {
                std::cout << RED << "Error: vertex " << u << " has " << in_same_partition
                          << " nodes in the same partition, and " << in_other_partition
                          << " in the other partition." << RESET << std::endl;
                return;
            }

            cut_size += in_other_partition;
        }
        assert(cut_size % 2 == 0);

        std::cout << GREEN << "Passed." << RESET << " Computed cut: " << cut_size / 2 << " edges. "
                  << std::endl;
    }
};

#endif /* ifndef MAX_CUT_HPP_ */
