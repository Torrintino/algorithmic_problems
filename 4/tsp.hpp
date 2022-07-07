#ifndef TSP_HPP_
#define TSP_HPP_

#include "problem.hpp"
#include "tsp-utils.hpp"

#include <random>
#include <vector>

//I ADDED THIS!!!
#include <algorithm>
//--------------

// Solution for TSP
using Cycle = std::vector<node>;

class TSP : public Problem<Graph, Cycle> {

    // To generate random nodes
    std::mt19937_64 gen;
    std::uniform_int_distribution<node> dist;

public:
    TSP(const Graph &G, double alpha = 100, double beta = 10, int seed = 42)
        : Problem(G, alpha, beta), dist(0, G.edges.size() - 1) {
        gen.seed(seed);
    }

    // Generates a new random Hamiltonian cycle.
    Cycle generate_initial_solution() override {
        
	    Cycle sol;
	    sol.reserve(inst.n);

	    for(node i = 0; i < inst.n; i++)
	    {
		sol.push_back(i);
	    }

	    std::shuffle(sol.begin(), sol.end(), gen);
	
	    return sol;
    }

    // Generates a new cycle by reverting a portion of the cycle.
    Cycle generate_random_neighbor(const Cycle &cycle) override {
		

	Cycle rnd_neighbor = cycle;

    	uint32_t x_i = dist(gen);
	uint32_t x_j = dist(gen);	
	
	if(x_i < x_j) std::reverse(rnd_neighbor.begin() + x_i, rnd_neighbor.begin() + x_j + 1 );
	else  std::reverse(rnd_neighbor.begin() + x_j, rnd_neighbor.begin() + x_i + 1 );
	
	return rnd_neighbor;
    }

    // Computes the cost of a cycle.
    double f(const Cycle &cycle) const override {
        
	    double cost = 0;

	    for(int i = 0; i < inst.n - 1; i++)
	    {
		    cost += inst[cycle.at(i)][cycle.at(i + 1)];
	    }

	    return cost;
    }

    // Returns true if the equilibrium for a fixed temperature is reached, false otherwise.
    bool equilibrium_reached(uint32_t iterations, uint32_t accepted_moves) const noexcept override {
        
	   if ( alpha * inst.n <= iterations ) return true;

	   if ( beta * inst.n <= accepted_moves ) return true;

	   return false;
    }

    // Returns true if the stopping condition is met, false otherwise.
    bool stopping_condition(uint32_t iterations, uint32_t accepted_moves) const noexcept override {
        
	    if ( (alpha * inst.n <= iterations) && accepted_moves == 0) return true;
	    else return false; 
    }
};

#endif // TSP_HPP_
