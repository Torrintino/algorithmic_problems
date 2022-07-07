#include "max-sat-utils.hpp"
#include "vnd.hpp"

int main(int argc, char **argv) {
    

    Formula f = SAT_IO::read_formula(20);
    int l = 45;
    VND_MAX_SAT vnd(f, l);
    std::vector<uint32_t> solutions = vnd.run_vnd(); 
	
    vector_to_file(solutions, "./results/formula20_l" + std::to_string(l));
    
    return 0;
}
