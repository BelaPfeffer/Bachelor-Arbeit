#include <sdsl/csa_bitcompressed.hpp>
#include <vector>
#include <string> 

void testCorrectness(const std::string& text, const std::string& kmer, const std::vector<uint64_t>& output);
    // Implement the testCorrectness function here
    // Compare the output with the expected output
    // If they match, print "Test passed"
    // Otherwise, print "Test f
void testRandomSequence(const std::string& text, int k);

std::string findRandSequence(const std::string& text, std::size_t length);

std::vector<std::string> findRandQueries (const std::string& text, std::size_t length, std::size_t num_sequences);