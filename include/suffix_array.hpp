#ifndef SUFFIXARRAY_H
#define SUFFIXARRAY_H
#include "results.hpp"
#include <string>
#include <vector>
#include <cstdint>          // for uint64_t
#include <unordered_map>
#include <sdsl/suffix_arrays.hpp>

class SuffixArray {
public:
    // Constructors / destructor
    explicit SuffixArray(const std::string& input);
    SuffixArray() = default;
    ~SuffixArray() = default;

    // Introspection / stats
    virtual size_t memoryUsageBytes(MemoryResults& mRes) const;

    // Accessors
    const std::vector<uint64_t>& getSuffixArray() const;
    const std::vector<int>&      getLCPArray() const;
    const std::string&           getText() const;

    // Debug printing
    void printSuffixArray() const;
    void printLCPArray() const;

    // Search interfaces
    std::vector<uint64_t> search(const std::string& pattern) const;
    std::pair<std::vector<uint64_t>, std::vector<uint64_t>>
    search_val_and_pos(const std::string& pattern) const;
    void save(const std::string& filename) const;
    static SuffixArray load(const std::string& filename);

protected:
    // Members
    std::string              text;
    std::vector<uint64_t>    suffixArray; // positions in text
    std::vector<int>         lcpArray;    // LCP between SA[i] and SA[i+1]

    // Builders (used by ctor)
    void buildSuffixArray();           // naive O(n log^2 n) version
    void buildSuffixArrayOptimized();  // O(n log n) prefix-doubling
    void buildLCPArray();              // Kasai's algorithm, O(n)
};

#endif // SUFFIXARRAY_H
