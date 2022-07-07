#ifndef VND_MAX_SAT_HPP_
#define VND_MAX_SAT_HPP_

#include <algorithm>
#include <vector>

// If needed, add further includes here

#include "max-sat-utils.hpp"

class VND_MAX_SAT final {
    const Formula &formula;
    const size_t l_max;
    Assignment best_solution;

    // Iterates over every subset of literals with the given size
         template <class Lambda>
	 void forEachSubsetOfLiterals(size_t size, Lambda lambda) {
     		 std::vector<uint32_t> index_vec;
	 	 index_vec.reserve(size);
		 
		 std::vector<bool> perm(formula.n);
		 for (size_t i = 0; i < size; ++i)
			 perm[i] = true;	 
		 do {
			 uint32_t i = 0;
			 do {
				 if (perm[i])
					 index_vec.push_back(i);
				 ++i;
			 } while (index_vec.size() < size);	 
			 lambda(index_vec);
			 index_vec.clear();
		 } while (std::prev_permutation(perm.begin(), perm.end()));
     	 }

    // Returns the best solution in the l-neighborhood of the given solution
    Assignment best_in_neighborhood(Assignment current_solution, uint32_t satisfied, size_t l) {
        
	uint32_t sat = satisfied;
	Assignment solution = current_solution;
	
	
	forEachSubsetOfLiterals(l, [&](const std::vector<uint32_t> &literals) {
				
			
			Assignment flip = solution;

			for(int i = 0; i<literals.size(); i++)
				flip[i] = !flip[i];

			uint32_t newSat = formula.satisfied_clauses(flip);

			if(newSat > sat)
			{
				solution = flip;
				sat = newSat;
			}

			});

	return solution;
	
    }

public:
    VND_MAX_SAT(const Formula &formula, size_t l_max = 5) : formula(formula), l_max(l_max) {}

    std::vector<uint32_t> run_vnd() {
        
	    
	    std::vector<uint32_t> solutions;

	    Assignment best_solution = formula.generate_random_assignment();
	    uint32_t satisfied = formula.satisfied_clauses(best_solution);

	    size_t l = 1;

	    while(l <= l_max)
	    {
		
		    Assignment newSolution = this->best_in_neighborhood(best_solution, satisfied, l);

		    uint32_t newSatisfied = formula.satisfied_clauses(newSolution);
		    solutions.push_back(newSatisfied);

		    if(newSatisfied > satisfied)
		    {
			    best_solution = newSolution;
			    satisfied = newSatisfied;
			    l=1;
		    }
		    else l++;

	    }


	    return solutions;
    }

    const Assignment &get_best_solution() const { return best_solution; }
};

#endif /* ifndef VND_MAX_SAT_HPP_ */
