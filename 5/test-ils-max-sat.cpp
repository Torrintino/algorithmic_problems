#include "max-sat-utils.hpp"
#include "max-sat.hpp"

int main(int argc, char **argv) {

	Formula f = SAT_IO::read_formula(150);
	Assignment a = f.generate_random_assignment();	
	MAX_SAT msat(f);

	std::vector<uint32_t> solution1 = msat.run_ils(a, MAX_SAT::Perturbation::FIXED, MAX_SAT::Convergence::ITERS );

	std::vector<uint32_t> solution2 = msat.run_ils(a, MAX_SAT::Perturbation::FIXED, MAX_SAT::Convergence::ADAPTIVE );

	std::vector<uint32_t> solution3 = msat.run_ils(a, MAX_SAT::Perturbation::DYNAMIC, MAX_SAT::Convergence::ITERS );

	std::vector<uint32_t> solution4 = msat.run_ils(a, MAX_SAT::Perturbation::DYNAMIC, MAX_SAT::Convergence::ADAPTIVE );

	std::string path1 ("./resultsILS/formula_150FI");
	std::string path2 ("./resultsILS/formula_150FA");
	std::string path3 ("./resultsILS/formula_150DI");
	std::string path4 ("./resultsILS/formula_150DA");
	
	vector_to_file(solution1, path1);
   	vector_to_file(solution2, path2);
	vector_to_file(solution3, path3);
	vector_to_file(solution4, path4);

       	return 0;
}
