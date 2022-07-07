#include <math.h>
#include "max-sat-utils.hpp"
#include "max-sat.hpp"

int main(int argc, char **argv) {
    int sizes[] = {
        50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150,
        160, 170, 180, 190, 200};
    std::vector<uint32_t> tries_vector = {1, 5, 10, 100, 1000};
    for(int tries : tries_vector) {
        std::vector<uint32_t> result;
        std::vector<uint32_t> expected_result;
        for(int len : sizes) {
            Formula f = SAT_IO::read_formula(len);
            uint32_t actual = 0;
            for(int i=0; i<tries; i++) {
                Assignment a = f.generate_random_assignment();
                actual += f.satisfied_clauses(a);
            }
            uint32_t expected = len * 10 / (1-pow(2, f.k));
            actual = actual / tries;
            result.push_back(actual);
            expected_result.push_back(expected);
        }
        vector_to_file(result, "actual" + std::to_string(tries) + ".bin");
        vector_to_file(expected_result,
                       "expected" + std::to_string(tries) + ".bin");
    }
    vector_to_file(tries_vector, "tries.bin");
    return 0;
}
