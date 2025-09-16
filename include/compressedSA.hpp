#pragma once

#include "hashValue.hpp"
#include "compressV2.hpp"
#include <unordered_map>
#include <vector>
#include <string>
#include <iostream>

class compressedSA
{   
private:
    std::unordered_map<uint64_t, hashValue> hashMap;
    std::vector<int> CSA;
    std::string text;

public:
    void printMap(uint64_t k)
{
    for (const auto& [key, value] : this -> hashMap) {

        std::cout << "Key: \"" << key << "\""
                  << ", Decoded_Key: " << decode_dna5(key,k) // Assuming k=2 for decoding
                  << ", cSAindex: " << value.cSAindex
                  << ", occurences: " << value.occurences
                  << ", lcp_interval_index: " << value.lcp_interval_index
                  << ", shift: " << value.shift
                  << ", refOcc: " << value.refOccurrences
                  << ", trace: " << value.traceback_key 
                  << ", processed: " << value.processed
                  << std::endl;
    }
}
    uint64_t encode_dna5(const std::string& kmer);
    std::string decode_dna5(uint64_t encoded, unsigned k);
    std::vector<int> findPattern(std::string& kmer, unsigned k);
    static compressedSA compute (const std::string& fastaData, const unsigned k);
    size_t memoryUsageBytes() const;



    compressedSA(std::unordered_map<uint64_t, hashValue> hashMap, std::vector<int> CSA, std::string text)
    {
        this -> hashMap = hashMap;
        this -> CSA = CSA;
        this -> text = text;
    }

    compressedSA (const std::string& fastaData, const unsigned k)
    {
       *this = compute(fastaData, k);
    }
    ~compressedSA() {}
};