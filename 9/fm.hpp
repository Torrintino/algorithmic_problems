#ifndef FM_HPP_
#define FM_HPP_

#include <algorithm>
#include <cassert>
#include <iostream>
#include <limits>
#include <random>
#include <vector>

#include <networkit/auxiliary/BucketPQ.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/graph/GraphTools.hpp>

// Include further headers if needed

using NetworKit::count;
using NetworKit::edgeweight;
using NetworKit::Graph;
using NetworKit::node;

// Generate a balanced random partition uniformly at random.
std::vector<bool> generate_random_partition(const Graph *g, int seed = 42) {
    std::mt19937_64 gen;
    gen.seed(seed);

    std::vector<bool> randomPartition(g->upperNodeIdBound());
    // Put ceil(n / 2) vertices in the first partition.
    std::fill(randomPartition.begin(),
              randomPartition.begin() + std::ceil(g->upperNodeIdBound() / 2.), true);
    std::shuffle(randomPartition.begin(), randomPartition.end(), gen);

    return randomPartition;
}

// Computes the cut for a given graph and for a given partition.
count compute_cut(const Graph *g, const std::vector<bool> &partition) {
    count cut = 0;
    g->forNodes([&](node u) {
        g->forNeighborsOf(u, [&](node v, edgeweight weight) {
            if (partition[u] != partition[v])
                cut += static_cast<count>(weight);
        });
    });
    return cut / 2;
}

class FM final {
    // Pointer to the input graph
    const Graph *g;
    std::vector<bool> partition; // Whether or not a node is in partition 1.
    std::vector<bool> blocked;   // Whether or not a node has already been swapped
    const int64_t max_degree;

    Aux::BucketPQ pq1, pq2;

    static constexpr int64_t inf = std::numeric_limits<int64_t>::max();

    // Returns -gain(u) (i.e., int_degree(u) - ext_degree(u)), so that on top of the PQ we have the
    // node with highest gain.
    int64_t key_of_node(node u) const {
        int64_t gain = 0;
        g->forNeighborsOf(u, [&](node v, edgeweight weight) {
            if (partition[u] != partition[v]) // External degree
                gain -= static_cast<int64_t>(weight);
            else // Internal degree
                gain += static_cast<int64_t>(weight);
        });

        return gain;
    }

    // All the vertices in V_1 are inserted into pq1, the vertices in V_2 into pq2.
    void fill_prio_queues() {
        g->forNodes([&](node u) { (partition[u] ? pq1 : pq2).insert(key_of_node(u), u); });
    }

    void clear_prio_queues() {
        while (pq1.size() > 0)
            pq1.extractMin();
        while (pq2.size() > 0)
            pq2.extractMin();
    }

    // Balances the current partition (if not balanced).
    void balance_partition() {
        
	    fill_prio_queues();
	    std::pair<int64_t,node> pair;

	    while(abs(pq1.size() - pq2.size()) > 1)
	    {
		    if(pq1.size() > pq2.size())
		    {
			    pair = pq1.extractMin();
			    partition[pair.second] = !partition[pair.second];
			    pq2.insert(key_of_node(pair.second), pair.second);
		    }
		    else
		    {
			    pair = pq2.extractMin();
			    partition[pair.second] = !partition[pair.second];
    			    pq1.insert(key_of_node(pair.second), pair.second);
		    }

		    g->forNeighborsOf(pair.second, [&](node v){
				    (partition[v] ? pq1 : pq2).changeKey(key_of_node(v),v);
				    });
	    }

	    clear_prio_queues();
	    
    }

public:
    // Takes an existing partition as input, balances it if needed.
    FM(const Graph *g, std::vector<bool> partition)
        : g(g), partition(std::move(partition)), blocked(g->upperNodeIdBound()),
          max_degree(NetworKit::GraphTools::maxWeightedDegree(*g)),
          pq1(g->upperNodeIdBound(), -max_degree, max_degree),
          pq2(g->upperNodeIdBound(), -max_degree, max_degree) {

        balance_partition();
    }

    // Generates a random partition.
    FM(const Graph *g, int seed = 42) : FM(g, generate_random_partition(g, seed)) {}

    // Returns the resulting partition.
    const std::vector<bool> &get_partition() const noexcept { return partition; }

    int64_t run_fm() {
	
	count initialCut = compute_cut(g, partition);
	count lastCut;
	count bestCut = initialCut;
	int bestCutIndex = 0;

	bool alternate;
	std::pair<int64_t,node> pair;
	
	int n = g->numberOfNodes();
	std::vector<node> move(n);
	

	int loopCount;

	//loop until no better solution found
	do
	{	
		loopCount = 0;

		blocked = std::vector<bool>(n);

		lastCut = bestCut;

		fill_prio_queues();
		
		
		while(pq1.size() > 0 || pq2.size() > 0)
		{	
			pair = (alternate ? pq1 : pq2).extractMin();
			alternate = !alternate;

			node u = pair.second;

			if(u >= inf) continue;

			partition[u] = !partition[u];
			move[loopCount] = u;
			blocked[u] = true;
			

			g->forNeighborsOf(u, [&](node v){
				if(blocked[v] == false) (partition[v] ? pq1 : pq2).changeKey(key_of_node(v), v);
				});

			count cut = compute_cut(g, partition);

			if(bestCut > cut)
			{
				bestCut = cut;
				bestCutIndex = loopCount;
			}		

			loopCount++;
		}

		if(bestCut < lastCut) for(int i = bestCutIndex + 1; i < n; i++) partition[move[i]] = !(partition[move[i]]);

		move.clear();
		clear_prio_queues();

	}while(bestCut < lastCut);

	return initialCut - lastCut;
    }
};

#endif /* ifndef FM_HPP_ */
