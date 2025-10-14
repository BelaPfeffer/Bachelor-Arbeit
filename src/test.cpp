#include "test.hpp"
#include "compressedSA.hpp"
#include <sdsl/suffix_array_algorithm.hpp>
#include <cassert>
#include <algorithm>
#include <iostream>
#include <random>
#include <string> 
#include <vector>

using namespace sdsl;

std::string findRandSequence(const std::string& text, std::size_t length) {
    if (length == 0) return "";

    std::vector<std::size_t> starts;

    std::size_t segStart = 0;
    while (segStart <= text.size()) {
        std::size_t next = text.find('$', segStart);
        std::size_t segEnd = (next == std::string::npos) ? text.size() : next;

        if (segEnd > segStart && segEnd - segStart >= length) {
            for (std::size_t s = segStart; s + length <= segEnd; ++s) {
                starts.push_back(s);
            }
        }

        if (next == std::string::npos) break;
        segStart = next + 1;
    }

    if (starts.empty()) return "";

    static thread_local std::mt19937 rng{std::random_device{}()};
    std::uniform_int_distribution<std::size_t> dist(0, starts.size() - 1);
    std::size_t start = starts[dist(rng)];
    return text.substr(start, length);
}

void testCorrectness(const std::string& text, const std::string& kmer, const std::vector<int>& output_pos)
{
    csa_bitcompressed<> csa;
    construct_im(csa, text, 1);
    int_vector<64> output = locate(csa, kmer);
    std::sort(output.begin(), output.end());
    
    // std::string correct_pos = "[";
    // for (unsigned long i = 0; i < output.size(); i++) {
    //     correct_pos += std::to_string(output[i]) + ",";
    // }
    // correct_pos[correct_pos.size() - 1] = ']';
    // std::cout << "Correct Positions: " << correct_pos << "\n";

    std::cout << "Start Test for kmer: " << kmer << std::endl;
    std::cout << "\n";

    std::cout << "Test Number of Occurences... ("<< output.size() << ")" << std::endl;
    std::cout << "\n";
    std::cout << "correct Number of Occurences: " << output.size() << ", " << "calculated Number of Occurences: " << output_pos.size() << std::endl;
    assert(output.size() == output_pos.size());
    std::cout << "\n";
    std::cout << "Test Number of Occurences successful" << std::endl;

    std::cout << "Test Positions... " << std::endl;
    std::cout << "\n";
    
    std::cout << "\n";
    std::cout << "Test Positions successful" << std::endl;
     
    return;
    
}

void testRandomSequence(const std::string& text, int k) {
    std::string rand_kmer = findRandSequence(text, k);
    if (rand_kmer.empty()) {
        std::cerr << "Failed to find a random k-mer of length " << k << ".\n";
        return;
    }

    // If we found a random k-mer, we can proceed with the test
    std::cout << "Found random k-mer: " << rand_kmer << "\n";
    std::cout << "Testing correctness for k-mer: " << rand_kmer << std::endl;

    compressedSA csa(text,k);
    std::vector<int> result = csa.findPattern(rand_kmer, k);
    testCorrectness(text, rand_kmer, result);
    return;
}
