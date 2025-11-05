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

std::vector<uint64_t> compressedSA::findPattern(std::string& kmer, unsigned k)
{
    if (kmer.size() != k) {
        throw std::invalid_argument("Pattern size must be equal to kmer size");
    }

    std::string retString = "[";

    std::vector<uint64_t> positions;
    bool isReference;
    bool isinCSA;

    uint64_t encoded_kmer = encode_dna5(kmer);
    hashValue curr_value = hashMap[encoded_kmer];
;

    isReference = (curr_value.refOccurrences != 0);
    isinCSA     = (curr_value.occurences != 0);

    // std::cout << "isReference: " << isReference << ", isinCSA: " << isinCSA << "\n";
    if (isinCSA) {
        unsigned occ       = curr_value.occurences;
        unsigned csa_index = curr_value.cSAindex;
        // std::cout << "csa_index: " << csa_index << ", occ: " << occ << "\n";
        for (unsigned long i = csa_index; i < csa_index + occ; ++i) {
            positions.emplace_back(CSA[i]);
            retString += std::to_string(CSA[i]) + ",";
        }
    }

    if (isReference) {
        unsigned trace          = curr_value.traceback_key;
        unsigned long refOcc    = curr_value.refOccurrences;
        int shift               = curr_value.shift;
        // std::cout << "trace: " << trace << ", refOcc: " << refOcc << ", shift: " << shift << "\n";

        for (unsigned long i = trace; i < trace + refOcc; ++i) {
            positions.emplace_back(CSA[i] + shift);
            retString += std::to_string(CSA[i] + shift) + ",";
        }
    }

    retString[retString.size() - 1] = ']';

    // std::cout << "Pattern kommt " << positions.size() << " mal vor in der Text, An Positionen: " << retString << "\n";

    return positions;
}



size_t compressedSA::memoryUsageBytes(MemoryResults& mRes) const {
    size_t totalMemory = 0;
    
    // Basis-Objektgröße (die Klasse selbst)
    totalMemory += sizeof(*this);
    
    size_t hashmapMEM = hashMap.size() * sizeof(std::pair<uint64_t, hashValue>);
    totalMemory += hashmapMEM;

    // Speicher für CSA vector
    size_t csaMEM = CSA.capacity() * sizeof(uint64_t);
    totalMemory += csaMEM;
    // Speicher für text string
    size_t textMEM = text.capacity() * sizeof(char);
    totalMemory += textMEM;
    
    std::cout << "Hashmap MEM (approx): " << hashmapMEM << " bytes\n";
    std::cout << "CSA size 8byte (approx): " << csaMEM << " bytes\n";
    std::cout << "CSA size 4byte (approx): " << csaMEM / 2 << " bytes\n";
    std::cout << "Text size (approx): " << textMEM << " bytes\n";

    mRes.setCSAData(csaMEM, hashmapMEM, CSA.size());

    return totalMemory;
}

compressedSA compressedSA::compute (const std::string& fastaData, const unsigned k) {
    std::unique_ptr<computeSA> csa = std::make_unique<computeSA>(fastaData, k);
    csa -> printAvgLCP();
    // csa -> printSuffixArray();
    // csa -> printIntervals(k);
    csa -> runCompression(k);
    // csa -> printMap(k);
    // csa -> printComputeSA(); 
    compressedSA e_csa = csa -> exportSA();
    // e_csa.printMap(k);

    return e_csa;
}

  void compressedSA::save(const std::string& filename) const {
        std::ofstream out(filename, std::ios::binary);
        if (!out) {
            throw std::runtime_error("Cannot open file for writing: " + filename);
        }
        
        // 1. Save text
        size_t text_size = text.size();
        out.write(reinterpret_cast<const char*>(&text_size), sizeof(text_size));
        out.write(text.data(), text_size);
        
        // 2. Save CSA vector
        size_t csa_size = CSA.size();
        out.write(reinterpret_cast<const char*>(&csa_size), sizeof(csa_size));
        out.write(reinterpret_cast<const char*>(CSA.data()), csa_size * sizeof(uint64_t));
        
        // 3. Save hashMap
        size_t map_size = hashMap.size();
        out.write(reinterpret_cast<const char*>(&map_size), sizeof(map_size));
        
        for (const auto& [key, value] : hashMap) {
            // Write key
            out.write(reinterpret_cast<const char*>(&key), sizeof(key));
            
            // Write hashValue (assuming it's POD - adjust if it has pointers/complex members)
            out.write(reinterpret_cast<const char*>(&value), sizeof(hashValue));
        }
        
        out.close();
        std::cout << "Saved compressedSA to " << filename << "\n";
    }
    
    // NEW: Load from disk
    compressedSA compressedSA::load(const std::string& filename) {
        std::ifstream in(filename, std::ios::binary);
        if (!in) {
            throw std::runtime_error("Cannot open file for reading: " + filename);
        }
        
        std::string text_loaded;
        std::vector<uint64_t> CSA_loaded;
        std::unordered_map<uint64_t, hashValue> hashMap_loaded;
        
        // 1. Load text
        size_t text_size;
        in.read(reinterpret_cast<char*>(&text_size), sizeof(text_size));
        text_loaded.resize(text_size);
        in.read(&text_loaded[0], text_size);
        
        // 2. Load CSA vector
        size_t csa_size;
        in.read(reinterpret_cast<char*>(&csa_size), sizeof(csa_size));
        CSA_loaded.resize(csa_size);
        in.read(reinterpret_cast<char*>(CSA_loaded.data()), csa_size * sizeof(uint64_t));
        
        // 3. Load hashMap
        size_t map_size;
        in.read(reinterpret_cast<char*>(&map_size), sizeof(map_size));
        
        for (size_t i = 0; i < map_size; ++i) {
            uint64_t key;
            hashValue value;
            
            in.read(reinterpret_cast<char*>(&key), sizeof(key));
            in.read(reinterpret_cast<char*>(&value), sizeof(hashValue));
            
            hashMap_loaded[key] = value;
        }
        
        in.close();
        
        std::cout << "Loaded compressedSA from " << filename << "\n";
        return compressedSA(hashMap_loaded, CSA_loaded, text_loaded);
    }