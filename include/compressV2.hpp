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
#include <sdsl/util.hpp>

using namespace sdsl;

struct lcp_interval
{
    typename csa_bitcompressed<>::size_type left; // Start index of the interval
    typename csa_bitcompressed<>::size_type right;   // End index of the interval
    typename csa_bitcompressed<>::size_type min_index; // Minimum LCP value in the interval
    unsigned priority;
    lcp_interval(typename csa_bitcompressed<>::size_type s, typename csa_bitcompressed<>::size_type e, typename csa_bitcompressed<>::size_type m) : left(s), right(e), min_index(m){};
    lcp_interval() : left(0), right(0), min_index(0), priority(0) {}; // Default constructor
};



struct hashValue
{   //=================================================================
    //KMER steht im komprimierten Suffix Array
    std::vector <unsigned long> cSAindex; // Index im komprimierten Suffix Array
    unsigned long occurences; // Anzahl der Vorkommen
    unsigned long lcp_interval_index; // Index des LCP Intervalls
    //=================================================================
    //REFERENZ auf anderes Pattern, von dem dieses abhängt
    int shift; // Verschiebung im Text (wenn -1 dann gibt keine Verschriebung d.h. nicht von anderem Pattern abhängig
    unsigned long refOccurrences;
    unsigned traceback_key; // Referenz auf das andere Pattern (kann nullptr sein, wenn es kein anderes Pattern gibt)
    bool processed;
               // sondern direkt in SA)

    
    hashValue()
    {
        cSAindex = std::vector <unsigned long>(); // noch nicht initialisiert
        occurences = 0; // noch nicht initialisiert
        shift = -1;  // noch nicht initialisiert
        refOccurrences = 0; // noch nicht 
        traceback_key = 0; // noch nicht initialisiert
        processed = false; // noch nicht initialisiert
    }

    void setValue (unsigned long cSAindex, unsigned long occurences)
    {
        this->cSAindex.push_back(cSAindex); 
        this->occurences += occurences;
        this -> processed = true;
    }

    void setReferenceValue(int shift, unsigned long refOccurrences, unsigned traceback_key, bool processed)
    {
        this->shift = shift;
        this->refOccurrences = refOccurrences;
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

    rank_support_v<> rankSupport; // Rank support for computeSuffix

    std::vector<std::optional<lcp_interval>> lcpIntervals; // Stores LCP intervals

    std::unordered_map<uint64_t, hashValue> hashMap;

    std::vector<int> compressedSA;

    void initComputeSuffix(unsigned k);

    void initLCPintervalsAndHashmap (const unsigned k);

    void initHash(const uint64_t& kmere, unsigned long lcp_interval_index);

    void readPattern (std::string& pattern, unsigned k);

    uint64_t encode_dna5(const std::string& kmer);

    std::string decode_dna5(uint64_t encoded, unsigned k);



public:
    void setBits(lcp_interval& iv, bool bitVal);

    void printMap(uint64_t k);

    bool filterSA (const int k, unsigned i);

    void printSuffixArray() const;

    lcp_bitcompressed<> get_lcp_array() const
    {
        return lcpArray;
    }

    size_t memoryUsageBytes() const;
    
    void printCompressedSA();

    hashValue getHashMapValue(const std::string& kmer);

    void setReferenceValue(const uint64_t encodedKmer, const int shift, unsigned long refOcurrences, const uint64_t traceback_key,  const bool processed);

    void setValue(const uint64_t encodedKmer, const unsigned long cSAindex, const unsigned long occurences);

    // void initHashMap(const unsigned k);

    void compression (const unsigned k, lcp_interval& interval);

   

    lcp_interval get_lcp_interval(unsigned i, unsigned k);

    void printIntervals(unsigned k);

    unsigned getInterval(unsigned i);

    unsigned calc_priority(lcp_interval& interval) const;

    void runCompression(const unsigned k);

    std::vector<int> findPattern(std::string& kmer, unsigned k);

    CompressedSA() = default;
    CompressedSA(const std::string& input, unsigned k) 
    {
        text = input;
        construct_im(suffixArray, text, 1);
        construct_im(lcpArray, text, 1);
        initComputeSuffix(k);
        util::init_support(rankSupport, &computeSuffix);
        initLCPintervalsAndHashmap(k);
    }

    
};