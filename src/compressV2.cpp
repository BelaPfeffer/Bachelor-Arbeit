#include "compressV2.hpp"
#include "suffix_array.hpp"
#include <iostream>
#include <algorithm>
#include <iterator>  

#include <sdsl/suffix_array_algorithm.hpp>

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
        if (value.processed == false) {
            continue; // Skip uninitialized entries
        }
        std::cout << "Key: \"" << key << "\""
                  << ", Decoded_Key: " << decode_dna5(key,k) // Assuming k=2 for decoding
                  << ", cSAindex: " << value.cSAindex.size()
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
        // std::cout << "Kmer: " << kmer << "\n";
        // std::cout << "pos: " << pos << "\n";

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

        unsigned long steps = i + 1;
        if(steps == suffixArray.size())
        {
            break;
        }
        // std::cout << "steps: " << steps << "\n";
        

        
        while (static_cast<unsigned>(lcpArray[steps]) >= k) {steps++; };
        // std::cout << "lcp loopend:" << std::endl;
   


        i = steps;

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

std::string CompressedSA::decode_dna5(uint64_t encoded, unsigned k)
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
    this -> hashMap[encodedKmer].setReferenceValue(shift, traceback_key, processed);
}


void CompressedSA::setValue(const uint64_t encodedKmer, const unsigned long cSAindex, const unsigned long occurences) // Index im komprimierten Suffix Array
{
   this ->hashMap[encodedKmer].setValue(cSAindex, occurences);
}

void CompressedSA::printSuffixArray() const 
{
        std::cout << "Suffix Array:\n";
        std::cout << "Index\tPosText\tLCP\tbitvektor\tSuffix\n";
        for (unsigned long i = 0; i < suffixArray.size(); i++) {
            std::cout << i << "\t"<<suffixArray[i] << "\t" << lcpArray[i] << "\t" << computeSuffix[i] << "\t" << text.substr(suffixArray[i]) << "\n";
        }
}




lcp_interval CompressedSA::get_lcp_interval(unsigned  i, unsigned k) 
{
    const sdsl::lcp_bitcompressed<>& lcp = lcpArray;
    unsigned long n = lcp.size() + 1; // Because LCP array has size n - 1
    unsigned long left = i;
    unsigned long right = i;
    unsigned long min_index = i;

    // Expand to the left
    while (left > 0 && lcp[left - 1] >= k) {
         if (lcp[left] < lcp[min_index]) {
            min_index = left;
        }
        --left;
    }

    // Expand to the right
    while (right < n - 1 && lcp[right] >= k) {
        if (lcp[right] < lcp[min_index]) {
            min_index = right;
        }
        ++right;
    }

    return lcp_interval(left, right, min_index);
}


//wird true wenn wert valid
bool CompressedSA::filterSA(const int k, unsigned i)
{
    unsigned long suffix = suffixArray[i];
    if (suffix + k > text.size()) {
        return false; // Suffix is too short to contain a k-mer of length k
    }
    else if( text.substr(suffix, k).find('$') != std::string::npos) {
        return false; // Suffix contains the end character '$'
    }
    return true;
}

void CompressedSA::initComputeSuffix(unsigned k)
{
    computeSuffix = bit_vector(suffixArray.size());
    for (unsigned long i = 0; i < suffixArray.size(); i++) {
        if (filterSA(k, i)) {
            computeSuffix[i] = 1; // Mark this suffix as valid
        } else {
            computeSuffix[i] = 0; // Mark this suffix as invalid
        }
    }
}

void CompressedSA::compression(const unsigned k, lcp_interval& interval)
{
    unsigned long current_csa_index = compressedSA.empty() ? 0 : compressedSA.size();
    unsigned pat_pos_index = suffixArray[interval.min_index]; // Index des Musters im text 
    unsigned long occurences = interval.right - interval.left + 1; // Anzahl der Vorkommen des Musters
    uint64_t kmer_old = encode_dna5(text.substr(pat_pos_index, k));
    setValue(kmer_old, current_csa_index, occurences); // das könnte probleme geben wenn compressed SA noch leer ist


    // ab hier nochmal gucken

    int shift = 1;

    unsigned text_index = pat_pos_index + shift + k - 1;

    for (unsigned long i = pat_pos_index + shift; text[text_index] != '$'; i ++) {

        lcp_interval temp_interval = lcp_interval();
        auto count = backward_search(suffixArray,0 , suffixArray.size() - 1, 
                                    text.begin() + i, text.begin() + i + k,
                                    temp_interval.left, temp_interval.right);

        unsigned kmer_new = encode_dna5(text.substr(i,k));

        while (count > occurences)
        {
            unsigned long x = 0;
            unsigned long y = 0;
            auto text_pos_kmer_new = suffixArray[temp_interval.left + x];
            auto text_pos_kmer_old = suffixArray[pat_pos_index + y];

            if (text_pos_kmer_new != text_pos_kmer_old)
            {
                setValue(kmer_new, compressedSA.size(), 1);
                y++;
                count--;
                continue;
            }
            x++;
            y++;
        }

        
        
        setReferenceValue(kmer_new, shift, hashMap[kmer_old].cSAindex[0], true);
        text_index ++;
        shift++;
    
    }

}
