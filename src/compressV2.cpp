#include "compressV2.hpp"
#include "suffix_array.hpp"
#include <iostream>
#include <algorithm>
#include <iterator>  

void CompressedSA::printCompressedSA()
{
    std::cout << "Compressed Suffix Array:\n";
    for (unsigned long i = 0; i < compressedSA.size(); i++) {
        std::cout << i << "\t" << compressedSA[i] << "\t" << this -> text.substr(compressedSA[i]) << "\n";
    }
    std::cout << "\n";
}

 size_t CompressedSA::memoryUsageBytes() const {
    size_t size = sizeof(*this); // Speicher für CompressedSA selbst

    // Dynamisch belegter Speicher im Vektor
    size += compressedSA.capacity() * sizeof(int);

    // Dynamisch belegter Speicher in der Hashmap
    size += hashMap.size() * (sizeof(uint64_t) + sizeof(hashValue));

    // Optional: Wenn load factor hoch ist, ist capacity > size:
    size += (hashMap.bucket_count() - hashMap.size()) * (sizeof(void*) + sizeof(size_t)); // ungefähr

    // Speicher aus der Basisklasse (SuffixArray)
    size += SuffixArray::memoryUsageBytes(); 

    return size;
}


    void CompressedSA::printMap(uint64_t k)
{
    for (const auto& [key, value] : this -> hashMap) {
        std::cout << "Key: \"" << key << "\""
                  << ", Decoded_Key: " << decode_dna5(key,k) // Assuming k=2 for decoding
                  << ", kSAindex: " << value.kSAindex
                  << ", occurences: " << value.occurences
                  << ", shift: " << value.shift
                  << ", trace: " << value.traceback_key 
                  << ", processed: " << value.processed
                  << std::endl;
    }
}

void CompressedSA::initHash(const uint64_t& kmer)
{
    this ->hashMap[kmer] = hashValue();   
}

void CompressedSA::initHashMap(const unsigned k)
{
    const std::vector<int>& suffixArray = this->suffixArray;
    const std::vector<int>& lcpArray = this->lcpArray;
    const std::string& text = this->text;
    unsigned long i = 0;

    while(i < suffixArray.size())
    {
        unsigned long pos = suffixArray[i];

        if (pos + k > text.size() )
        {
            // std::cout << "Kmer is shorter than k: " << kmer_old << "\n";
            i++;
            continue;
        }

         std::string kmer = text.substr(pos, k);

        if(kmer.find('$') != std::string::npos)
        {
            i++;
            continue;
        }

        uint64_t encodedKmer = encode_dna5(kmer);
        this -> initHash(encodedKmer);
        for (unsigned long j = 0; j < kmer.size(); j++)
        {
            std::cout << kmer[j];
        }
        std::cout << "\n";

        unsigned long steps = i;
        

        //jumps to next Bucket
        while (static_cast<unsigned>(lcpArray[steps]) >= k) steps++;
        
        
        i = steps + 1;
        
    }
    std::cout << "Map size: " << hashMap.size() << "\n";
}


uint64_t CompressedSA::encode_dna5(const std::string& kmer)
{
    uint64_t encoded = 0;
    for (char c : kmer) {
        encoded <<= 3; // Shift left by  bits
        switch (c) {
            case 'A': encoded |= 0; break; // 000
            case 'C': encoded |= 1; break; // 001
            case 'G': encoded |= 2; break; // 010
            case 'T': encoded |= 3; break; // 011
            case '$': encoded |= 4; break; // 100 (special end character)
            default: throw std::invalid_argument("Invalid character in DNA sequence");
        }
    }
    return encoded;
}

std::string CompressedSA::decode_dna5(uint64_t encoded, int k)
{
    std::string decoded;
    for (int i = 0; i < k; ++i) {
        int base = encoded & 0x07; // Get the last 3 bits
        switch (base) {
            case 0: decoded += 'A'; break;
            case 1: decoded += 'C'; break;
            case 2: decoded += 'G'; break;
            case 3: decoded += 'T'; break;
            case 4: decoded += '$'; break; // Special end character
            default: throw std::invalid_argument("Invalid base in encoded DNA sequence");
        }
        encoded >>= 3; // Shift right by 3 bits
    }
    std::reverse(decoded.begin(), decoded.end()); // Reverse to get the correct order
    return decoded;
}

hashValue CompressedSA::getHashMapValue(const std::string& kmer)
{
    uint64_t encodedKmer = encode_dna5(kmer);
    return this->hashMap[encodedKmer];
}

void CompressedSA::setReferenceValue(const uint64_t encodedKmer, const int shift, const uint64_t traceback_key, const bool processed)
{
    this->hashMap[encodedKmer].traceback_key = traceback_key;
    
    this->hashMap[encodedKmer].shift = shift;

    this->hashMap[encodedKmer].processed = processed;
}


void CompressedSA::setValue(const uint64_t encodedKmer, const unsigned long kSAindex, const unsigned long occurences) // Index im komprimierten Suffix Array
{
    this->hashMap[encodedKmer].kSAindex = kSAindex;

    this->hashMap[encodedKmer].occurences = occurences;
}


