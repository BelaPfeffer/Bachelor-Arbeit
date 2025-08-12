#include "test.hpp"
#include <sdsl/suffix_array_algorithm.hpp>
#include <cassert>
#include <iostream>

using namespace sdsl;

void testCorrectness(const std::string& text, const std::string& kmer, const std::vector<int>& output_pos)
{
    csa_bitcompressed<> csa;
    construct_im(csa, text, 1);
    int_vector<64> output = locate(csa, kmer);

    std::cout << "Start Test for kmer: " << kmer << std::endl;
    std::cout << "\n";

    std::cout << "Test Number of Occurences... ("<< output.size() << ")" << std::endl;
    std::cout << "\n";
    assert(output.size() == output_pos.size());
    std::cout << "\n";
    std::cout << "Test Number of Occurences successful" << std::endl;

    std::cout << "Test Positions... " << std::endl;
    std::cout << "\n";
    for (unsigned i = 0; i < output.size(); ++i)
    {
        std::cout << "correct Position: " << output[i] << ", " << "calculated Position: " << output_pos[i] << std::endl;
        assert(output[i] == output_pos[i]);
    }
    std::cout << "\n";
    std::cout << "Test Positions successful" << std::endl;
     
    return;
    
}