#include <map>
#include <tuple>
#include <string>
#include "suffix_array.hpp"

struct hashValue
{
    unsigned long kSAindex; // Index im komprimierten Suffix Array
    unsigned long occurences; // Anzahl der Vorkommen
    int shift; // Verschiebung im Text (wenn -1 dann gibt keine Verschriebung d.h. nicht von anderem Pattern abhängig
    hashValue* traceback;
    bool processed;
               // sondern direkt in SA)

    
    hashValue()
    {
        kSAindex = 0; // noch nicht initialisiert
        occurences = 0; // noch nicht initialisiert
        shift = -1;  // noch nicht initialisiert
        traceback = nullptr; // noch nicht initialisiert
        processed = false; // noch nicht initialisiert
    }

    
};





class CompressedSA : public SuffixArray
{
private:
    std::unordered_map<std::string, hashValue> hashMap;
    std::vector<int> compressedSA;

    std::pair<unsigned long, int> findLCPmaxmin(const int k);

    void readPattern (std::string& pattern, int k);



public:
    void printMap();

    void printCompressedSA();

    void initHashMap(const int k);

    void compression (const int k);

    void initHash(const std::string& kmere);

    CompressedSA() = default;
    CompressedSA(const std::string& input) : SuffixArray(input) {}
    ~CompressedSA() = default;
};