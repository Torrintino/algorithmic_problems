#ifndef MAX_SAT_HPP_
#define MAX_SAT_HPP_

#include "max-sat-utils.hpp"

#include <cassert>
#include <limits>
#include <random>
#include <unordered_set>
#include <vector>

constexpr auto none = std::numeric_limits<uint32_t>::max();

class MAX_SAT {
    // Instance of max-sat.
    const Formula &formula;

    std::mt19937_64 gen;

    // Parameters of the algorithm. In your analysis try to change some of them
    // and see how they affect the results.
    const uint32_t max_iters_ils = 30, max_iters_with_no_new_best = 10;
    const double fixed_perturb_p = 0.1, dyn_perturb_p_start = 0.5, dyn_perturb_p_mult = 0.9;
    const uint32_t max_gsat_iters = 100;

    // Runs gsat on the current solution
    void run_gsat(Assignment &assignment) {
        uint32_t satisfied_clauses = formula.satisfied_clauses(assignment);

        const auto literal_to_flip = [&]() -> uint32_t {
            std::vector<uint32_t> candidates;
            for (uint32_t i = 0; i < formula.n; ++i) {
                assignment[i] = !assignment[i];
                uint32_t newly_satisfied_clauses = formula.satisfied_clauses(assignment);
                assignment[i] = !assignment[i];
                if (newly_satisfied_clauses > satisfied_clauses) {
                    candidates = {i};
                    satisfied_clauses = newly_satisfied_clauses;
                } else if (!candidates.empty() && newly_satisfied_clauses == satisfied_clauses)
                    candidates.push_back(i);
            }

            if (candidates.empty())
                return none;

            const size_t random_candidate =
                std::uniform_int_distribution<size_t>{0, candidates.size() - 1}(gen);
            return candidates[random_candidate];
        };

        uint32_t iters_done = 0;

        bool in_local_optimum;
        do {
            if (satisfied_clauses == formula.clauses.size())
                return; // We are in a global optimum

            in_local_optimum = true;
            const uint32_t to_flip = literal_to_flip();
            if (to_flip != none) {
                assignment[to_flip] = !assignment[to_flip];
                in_local_optimum = false;
            }
        } while (!in_local_optimum && (max_gsat_iters == 0 || ++iters_done < max_gsat_iters));
    }

    // Flips ceil(n * p) random literals in the given assignment, where n is the number of literals.
    void flip_random_literals(Assignment &assignment, double p) {
        
	    double perturb = ceil(assignment.size()*p);

	    for(int i = 0; i < perturb; i++)
	    {
		    size_t c = std::uniform_int_distribution<size_t>{0, assignment.size() - 1}(gen);

		    assignment[c] = !assignment[c];
	    }
    }

public:
    MAX_SAT(const Formula &formula, uint32_t seed = 42) : formula(formula) { gen.seed(seed); }

    // Perturbation strategies
    enum class Perturbation { FIXED, DYNAMIC };

    // Convergence criteria
    enum class Convergence { ITERS, ADAPTIVE };

    // Runs ILS on the given assignment; follows the given perturbation strategy
    // and convergence criterion.
    std::vector<uint32_t> run_ils(Assignment &assignment,
                                  Perturbation perturbation = Perturbation::FIXED,
                                  Convergence convergence = Convergence::ITERS) {
        
	    Assignment a = assignment; 

	    uint32_t iterCount = 0;
	    uint32_t bestSolution = formula.satisfied_clauses(a);

	    std::vector<uint32_t> solutions = {};

	    double p;
	    if(perturbation == Perturbation::FIXED) p = fixed_perturb_p;
	    else p = dyn_perturb_p_start;

	    while( (convergence == Convergence::ITERS && iterCount <= max_iters_ils) || (convergence == Convergence::ADAPTIVE && iterCount <= max_iters_with_no_new_best) )
	    {

		flip_random_literals(a, p);

		run_gsat(a);
		
		uint32_t solution = formula.satisfied_clauses(a);
		solutions.push_back(solution);
		if(solution > bestSolution)
		{
			bestSolution = solution;
			if(convergence == Convergence::ADAPTIVE)
			{	
				iterCount = 0;	
			}
		}


		//Update p
		if(perturbation == Perturbation::DYNAMIC) p *= dyn_perturb_p_mult;

		iterCount++;
	    }
	   
	    
	    return solutions;

    }
};

#endif /* ifndef MAX_SAT_HPP_ */
