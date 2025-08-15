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