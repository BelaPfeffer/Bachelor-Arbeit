#include <sdsl/csa_bitcompressed.hpp>
#include <vector>
#include <string> 

void testCorrectness(const std::string& text, const std::string& kmer, const std::vector<int>& output);
    // Implement the testCorrectness function here
    // Compare the output with the expected output
    // If they match, print "Test passed"
    // Otherwise, print "Test f