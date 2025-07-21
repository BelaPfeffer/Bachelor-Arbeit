#include <map>
#include <tuple>
#include <string>
#include "suffix_array.hpp"

struct hashValue
{   //=================================================================
    //KMER steht im komprimierten Suffix Array
    unsigned long kSAindex; // Index im komprimierten Suffix Array
    unsigned long occurences; // Anzahl der Vorkommen
    //=================================================================
    //REFERENZ auf anderes Pattern, von dem dieses abhängt
    int shift; // Verschiebung im Text (wenn -1 dann gibt keine Verschriebung d.h. nicht von anderem Pattern abhängig
    uint64_t traceback_key; // Referenz auf das andere Pattern (kann nullptr sein, wenn es kein anderes Pattern gibt)
    bool processed;
               // sondern direkt in SA)

    
    hashValue()
    {
        kSAindex = 0; // noch nicht initialisiert
        occurences = 0; // noch nicht initialisiert
        shift = -1;  // noch nicht initialisiert
        traceback_key = 0; // noch nicht initialisiert
        processed = false; // noch nicht initialisiert
    }

    
};





class CompressedSA : public SuffixArray
{
private:
    std::unordered_map<uint64_t, hashValue> hashMap;

    std::vector<int> compressedSA;

    void readPattern (std::string& pattern, int k);

    uint64_t encode_dna5(const std::string& kmer);

    std::string decode_dna5(uint64_t encoded, int k);



public:
    void printMap();

    void printCompressedSA();

    hashValue getHashMapValue(const std::string& kmer);

    void setReferenceValue(const uint64_t encodedKmer, const int shift, const uint64_t traceback_key,  const bool processed);

    void setValue(const uint64_t encodedKmer, const unsigned long kSAindex, const unsigned long occurences);

    void initHashMap(const unsigned k);

    void compression (const int k);

    void initHash(const uint64_t& kmere);

    CompressedSA() = default;
    CompressedSA(const std::string& input) : SuffixArray(input) {}
    ~CompressedSA() = default;
};