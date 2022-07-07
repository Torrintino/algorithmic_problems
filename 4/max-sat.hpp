#ifndef MAX_SAT_HPP_
#define MAX_SAT_HPP_

#include "max-sat-utils.hpp"
#include "problem.hpp"

#include <random>
#include <vector>

class MAX_SAT : public Problem<Formula, Assignment> {

    // To pick random literals.
    std::mt19937_64 gen;
    std::uniform_int_distribution<uint32_t> dist;

public:
    MAX_SAT(const Formula &formula, double alpha = 100, double beta = 10, int seed = 42)
        : Problem(formula, alpha, beta), dist(0, formula.n - 1) {
        gen.seed(seed);
    }

    // Generates a new random assignment of the literals.
    Assignment generate_initial_solution() override {
       
	   Assignment initial;
	   initial.reserve(inst.n);

	   int prob50 = inst.n/2;

	   for(int i = 0; i < inst.n; i++)
	   {		
		   if(dist(gen) > prob50) initial.push_back(true);
		   else initial.push_back(false);
	   }

	   return initial;
	    
    };

    // Generates a new assignment by flipping a literal of the given assignment.
    Assignment generate_random_neighbor(const Assignment &a) override {
        
	    Assignment neighbor = a;

	    int rnd = dist(gen);

	    neighbor.at(rnd) = !neighbor.at(rnd);

	    return neighbor;
    }

    // Cost of the current assignment i.e., the number of UNsatisfied clauses in the formula.
    double f(const Assignment &a) const override {
        
	    double unsat = 0;


	    for(Clause clause : inst.clauses)
	    {
		    bool satisfied = false;

		    for(Literal x : clause)
		    {

			    if(a.at(x.i) == !x.neg )
			    {
				    satisfied = true;
				    break;
			    }
		    }

		    if(satisfied == false) unsat++;
	    }

	    return unsat;
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

#endif /* ifndef MAX_SAT_HPP_ */
