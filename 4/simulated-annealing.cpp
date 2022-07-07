#ifndef SIMULATED_ANNEALING_CPP_
#define SIMULATED_ANNEALING_CPP_

#include <cmath>
#include <random>

template <class ProblemType, class Instance, class Solution>
class SimulatedAnnealing {

    // Instance of the problem to be solved
    ProblemType problem;

    // Best solution found
    Solution best_solution;

    // Current temperature and parameter for the geometric cool down strategy
    double temperature, cool_down_param;

    // To determine the acceptance of a worse solution
    std::mt19937_64 gen;
    std::uniform_real_distribution<double> dist;

    // Returns true if the new solution (f_j) should be accepted, false otherwise
    bool accept(double f_i, double f_j) {
        
	    if( f_j <= f_i ) return true;
	    

	    double prob = exp(-((f_j - f_i)/temperature));
	    if (dist(gen) < prob) return true;

	    return false;
    }

    // Updates the temperature using the geometric cool down strategy
    void cool_down() noexcept {
        
	    temperature = temperature * cool_down_param;
    }

public:
    SimulatedAnnealing(const Instance &inst, double initial_temperature,
                       double cool_down_param = 0.9, int seed = 42)
        : problem(inst), temperature(initial_temperature), cool_down_param(cool_down_param),
          dist(0, 1) {
        gen.seed(seed);
    }

    // Runs the simulated annealing algorithm
    void run() {
        
	
	    uint32_t count_inner = 0;
	    uint32_t count_acc = 0;

	    best_solution = problem.generate_initial_solution();

	    while(!problem.stopping_condition(count_inner, count_acc) && temperature >= 0.001 )
	    {

		    count_inner = 0;
		    count_acc = 0;
		
		    while(!problem.equilibrium_reached(count_inner, count_acc))
		    {

			    Solution new_solution = problem.generate_random_neighbor(best_solution);

			    if(accept(problem.f(best_solution), problem.f(new_solution)))
			    {
				    best_solution = new_solution;
				    count_acc++;
			    }

			    count_inner++; 

		    }
	
		    cool_down();
	    }
	
    }

    // Returns a const ref to the best solution found
    const Solution &get_solution() const {
	    return best_solution;
    }
};

#endif
