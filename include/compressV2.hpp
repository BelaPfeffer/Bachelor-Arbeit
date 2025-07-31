#include <map>
#include <tuple>
#include <string>
#include <utility> 
#include "suffix_array.hpp"

#include <sdsl/csa_bitcompressed.hpp>
#include <sdsl/suffix_array_algorithm.hpp>
#include <sdsl/lcp_bitcompressed.hpp>
#include <sdsl/construct_lcp.hpp>
#include <sdsl/bit_vectors.hpp>

using namespace sdsl;

struct lcp_interval
{
    unsigned long left; // Start index of the interval
    unsigned long right;   // End index of the interval
    unsigned long min_index; // Minimum LCP value in the interval
    lcp_interval(unsigned long s, unsigned long e, unsigned long m) : left(s), right(e), min_index(m){};
};



struct hashValue
{   //=================================================================
    //KMER steht im komprimierten Suffix Array
    unsigned long cSAindex; // Index im komprimierten Suffix Array
    unsigned long occurences; // Anzahl der Vorkommen
    //=================================================================
    //REFERENZ auf anderes Pattern, von dem dieses abhängt
    int shift; // Verschiebung im Text (wenn -1 dann gibt keine Verschriebung d.h. nicht von anderem Pattern abhängig
    unsigned traceback_key; // Referenz auf das andere Pattern (kann nullptr sein, wenn es kein anderes Pattern gibt)
    bool processed;
               // sondern direkt in SA)

    
    hashValue()
    {
        cSAindex = 0; // noch nicht initialisiert
        occurences = 0; // noch nicht initialisiert
        shift = -1;  // noch nicht initialisiert
        traceback_key = 0; // noch nicht initialisiert
        processed = false; // noch nicht initialisiert
    }

    void setValue (unsigned long cSAindex, unsigned long occurences)
    {
        this->cSAindex = cSAindex;
        this->occurences = occurences;
    }

    void setReferenceValue(int shift, unsigned traceback_key, bool processed)
    {
        this->shift = shift;
        this->traceback_key = traceback_key;
        this->processed = processed;
    }

    
};





class CompressedSA : public SuffixArray
{
private:
    std::string text; // Originaltext

    csa_bitcompressed<> suffixArray;

    lcp_bitcompressed<> lcpArray;

    bit_vector computeSuffix;

    std::unordered_map<uint64_t, hashValue> hashMap;

    std::vector<int> compressedSA;

    void initComputeSuffix(unsigned k);

    void readPattern (std::string& pattern, unsigned k);

    uint64_t encode_dna5(const std::string& kmer);

    std::string decode_dna5(uint64_t encoded, unsigned k);



public:
    void printMap(uint64_t k);

    bool filterSA (const int k, unsigned i);

    void printSuffixArray() const;

    size_t memoryUsageBytes() const;
    
    void printCompressedSA();

    hashValue getHashMapValue(const std::string& kmer);

    void setReferenceValue(const uint64_t encodedKmer, const int shift, const uint64_t traceback_key,  const bool processed);

    void setValue(const uint64_t encodedKmer, const unsigned long cSAindex, const unsigned long occurences);

    void initHashMap(const unsigned k);

    void compression (const unsigned k, lcp_interval& interval);

    void initHash(const uint64_t& kmere);

    lcp_interval get_lcp_interval(unsigned i, unsigned k);

    CompressedSA() = default;
    CompressedSA(const std::string& input, unsigned k) 
    {
        text = input;
        construct_im(suffixArray, text, 1);
        construct_im(lcpArray, text, 1);
        initComputeSuffix(k);
        initHashMap(k);
        lcp_interval interval = get_lcp_interval( 0, k);

    }
};