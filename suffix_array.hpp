#ifndef SUFFIXARRAY_H
#define SUFFIXARRAY_H

#include <string>
#include <vector>

class SuffixArray {
public:
    // Constructor: Takes the input string
    SuffixArray(const std::string& input);
    SuffixArray() = default;
    ~SuffixArray() = default;

    // Public functions to get results or print them
    const std::vector<int>& getSuffixArray() const;
    const std::vector<int>& getLCPArray() const;
    const std::string& getText() const ;
    void printSuffixArray() const;
    void printLCPArray() const;
    
    // Public search function
    std::vector<int> search(const std::string& pattern) const;

protected:
    // Member variables
    std::string text;
    std::vector<int> suffixArray;
    std::vector<int> lcpArray;

    // Private helper functions used by the constructor
    void buildSuffixArray();           // The simple O(n^2 log n) version
    void buildSuffixArrayOptimized();  // The O(n log n) version
    void buildLCPArray();              // Kasai's algorithm
};

#endif // SUFFIXARRAY_H