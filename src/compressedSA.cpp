#include "compressedSA.hpp"
#include <string>



uint64_t compressedSA::encode_dna5(const std::string& kmer)
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

std::string compressedSA::decode_dna5(uint64_t encoded, unsigned k)
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

std::vector<int> compressedSA::findPattern(std::string& kmer, unsigned k)
{
    if(kmer.size() != k) throw std::invalid_argument("Pattern size must be equal to kmer size");

    std::string retString = "[";

    std::vector<int> positions;
    bool isReference;
    bool isinCSA;

    int encoded_kmer = encode_dna5(kmer);
    hashValue curr_value = hashMap[encoded_kmer];

    isReference = (curr_value.refOccurrences != 0);
    isinCSA = (curr_value.occurences != 0);
    std::cout << "isReference: " << isReference << ", isinCSA: " << isinCSA << "\n";
    if(isinCSA)
    {   
        unsigned occ = curr_value.occurences;
        unsigned csa_index = curr_value.cSAindex; 
        std::cout << "csa_index: " << csa_index << ", occ: " << occ << "\n";
        for (unsigned long i = csa_index; i < csa_index + occ; i++)
        {
            positions.emplace_back(CSA[i]);
            retString += std::to_string(CSA[i]) + ",";
        }
    }

    if (isReference)
    {
        unsigned trace = curr_value.traceback_key;
        unsigned long refOcc = curr_value.refOccurrences;
        int shift = curr_value.shift;
        // std::cout << "trace: " << trace << ", refOcc: " << refOcc << ", shift: " << shift << "\n";

        for (unsigned long i = trace; i < trace + refOcc; i++)
        {
            positions.emplace_back(CSA[i] + shift);
            retString += std::to_string(CSA[i] + shift) + ",";
        }
    }
    
    retString[retString.size() - 1] = ']';

    std::cout << "Pattern kommt " << positions.size() << " mal vor in der Text, An Positionen: " << retString << "\n";
    
    return positions;
}

size_t compressedSA::memoryUsageBytes() const {
    size_t totalMemory = 0;
    
    // Basis-Objektgröße (die Klasse selbst)
    totalMemory += sizeof(*this);
    
    totalMemory += hashMap.size() * sizeof(std::pair<uint64_t, hashValue>);
    
    // Speicher für CSA vector
    totalMemory += CSA.capacity() * sizeof(int);
    
    // Speicher für text string
    totalMemory += text.capacity() * sizeof(char);
    
    return totalMemory;
}

compressedSA compressedSA::compute (const std::string& fastaData, const unsigned k) {
    std::unique_ptr<computeSA> csa = std::make_unique<computeSA>(fastaData, k);
    // csa -> printSuffixArray();
    // csa -> printIntervals(k);
    // csa -> runCompression(k);
    // csa -> printMap(k);
    // csa -> printComputeSA();
    compressedSA e_csa = csa -> exportSA();
    return e_csa;
}