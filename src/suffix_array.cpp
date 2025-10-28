// suffix_array.cpp
#include "suffix_array.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

// --- Naive O(n log^2 n) build (kept for completeness) ---
void SuffixArray::buildSuffixArray() {
    const size_t n = text.size();
    suffixArray.resize(n);

    for (size_t i = 0; i < n; ++i)
        suffixArray[i] = static_cast<uint64_t>(i);

    std::sort(suffixArray.begin(), suffixArray.end(),
              [this](uint64_t a, uint64_t b) {
                  return text.substr(static_cast<size_t>(a)) <
                         text.substr(static_cast<size_t>(b));
              });
}

// Approx memory usage of the object and its dynamic buffers
size_t SuffixArray::memoryUsageBytes() const {
    size_t totalMemory = 0;

    // Base object size (the class itself)
    totalMemory += sizeof(*this);

    // Text buffer (uses capacity, like your compressedSA)
    const size_t textMEM = text.capacity() * sizeof(char);
    totalMemory += textMEM;

    // Suffix Array vector
    const size_t saMEM = suffixArray.capacity() * sizeof(uint64_t);
    totalMemory += saMEM;

    // LCP Array vector (kept as 'int' to match the header)
    const size_t lcpMEM = lcpArray.capacity() * sizeof(int);
    totalMemory += lcpMEM;

    // Verbose breakdown (approximate, ignores allocator/metadata overhead)
    std::cout << "Object size (approx): " << sizeof(*this) << " bytes\n";
    std::cout << "Text capacity (approx): " << textMEM << " bytes\n";
    std::cout << "SA capacity (approx): " << saMEM << " bytes\n";
    std::cout << "LCP capacity (approx): " << lcpMEM << " bytes\n";

    return totalMemory;
}


// --- O(n log n) prefix-doubling SA build (recommended) ---
void SuffixArray::buildSuffixArrayOptimized() {
    const size_t n = text.size();
    suffixArray.resize(n);

    // 32-bit ranks are sufficient and cache-friendly
    std::vector<int32_t> rank(n), tempRank(n);

    for (size_t i = 0; i < n; ++i) {
        rank[i] = static_cast<unsigned char>(text[i]); // avoid char sign issues
        suffixArray[i] = static_cast<uint64_t>(i);
    }

    for (size_t k = 1; k < n; k <<= 1) {
        std::sort(suffixArray.begin(), suffixArray.end(),
                  [&](uint64_t a, uint64_t b) {
                      if (rank[a] != rank[b]) return rank[a] < rank[b];
                      int32_t ra = (a + k < n) ? rank[a + k] : -1;
                      int32_t rb = (b + k < n) ? rank[b + k] : -1;
                      return ra < rb;
                  });

        tempRank[suffixArray[0]] = 0;
        for (size_t i = 1; i < n; ++i) {
            const uint64_t prev = suffixArray[i - 1];
            const uint64_t curr = suffixArray[i];
            const bool diff = (rank[prev] != rank[curr]) ||
                              ((prev + k < n ? rank[prev + k] : -1) !=
                               (curr + k < n ? rank[curr + k] : -1));
            tempRank[curr] = tempRank[prev] + (diff ? 1 : 0);
        }
        rank.swap(tempRank);

        // Early exit if all ranks are 0..n-1
        if (rank[suffixArray.back()] == static_cast<int32_t>(n - 1)) break;
    }
}

// --- Kasai LCP in O(n) ---
void SuffixArray::buildLCPArray() {
    const size_t n = text.size();
    if (n == 0) { lcpArray.clear(); return; }
    lcpArray.assign(n > 0 ? n - 1 : 0, 0);

    std::vector<size_t> invSA(n);
    for (size_t i = 0; i < n; ++i)
        invSA[static_cast<size_t>(suffixArray[i])] = i;

    size_t lcp = 0;
    for (size_t i = 0; i < n; ++i) {
        const size_t r = invSA[i];
        if (r == n - 1) { lcp = 0; continue; } // last SA entry has no next neighbor
        const size_t j = static_cast<size_t>(suffixArray[r + 1]);
        while (i + lcp < n && j + lcp < n && text[i + lcp] == text[j + lcp]) ++lcp;
        lcpArray[r] = static_cast<int>(lcp); // keep storage type 'int'
        if (lcp) --lcp;
    }
}

// --- Ctor ---
SuffixArray::SuffixArray(const std::string& input) : text(input) {
    if (!text.empty()) {
        buildSuffixArrayOptimized();
        buildLCPArray();
    }
}

// --- Getters ---
const std::vector<uint64_t>& SuffixArray::getSuffixArray() const { return suffixArray; }
const std::string& SuffixArray::getText() const { return text; }
const std::vector<int>& SuffixArray::getLCPArray() const { return lcpArray; }

// --- Debug printing (fixed bounds) ---
void SuffixArray::printSuffixArray() const {
    std::cout << "Suffix Array:\n";
    std::cout << "Index\tPosText\tLCP\tSuffix\n";
    const size_t n = suffixArray.size();
    for (size_t i = 0; i < n; ++i) {
        const int lcp = (i < lcpArray.size()) ? lcpArray[i] : 0; // safe for last row
        std::cout << i << "\t" << suffixArray[i] << "\t" << lcp
                  << "\t" << text.substr(static_cast<size_t>(suffixArray[i])) << "\n";
    }
}

void SuffixArray::printLCPArray() const {
    std::cout << "\nLCP Array:\n";
    std::cout << "Index\tLCP\tSuffix 1\t\tSuffix 2\n";
    if (suffixArray.size() < 2) return;
    for (size_t i = 0; i < lcpArray.size(); ++i) { // lcpArray has size n-1
        std::cout << i << "\t" << lcpArray[i] << "\t"
                  << text.substr(static_cast<size_t>(suffixArray[i])) << "\t\t"
                  << text.substr(static_cast<size_t>(suffixArray[i + 1])) << "\n";
    }
}

// --- FIXED Binary search for pattern matching ---
std::vector<uint64_t> SuffixArray::search(const std::string& pattern) const {
    std::vector<uint64_t> result;
    if (pattern.empty() || suffixArray.empty()) return result;

    const char* T = text.data();
    const size_t N = text.size();
    const size_t m = pattern.size();

    // Find the first suffix that is >= pattern (could start with pattern)
    auto first_it = std::lower_bound(
        suffixArray.begin(), suffixArray.end(), pattern,
        [&](uint64_t sa_pos, const std::string& p) {
            // Compare suffix at sa_pos with pattern
            size_t tlen = N - sa_pos;
            size_t cmp_len = std::min(m, tlen);
            int c = std::char_traits<char>::compare(T + sa_pos, p.data(), cmp_len);
            if (c != 0) return c < 0;
            // If equal up to cmp_len: suffix < pattern only if suffix is shorter than pattern
            return tlen < m;
        });

    // Find the first suffix that does NOT start with pattern
    // We compare only the first m characters
    auto last_it = first_it;
    while (last_it != suffixArray.end()) {
        uint64_t sa_pos = *last_it;
        size_t tlen = N - sa_pos;
        size_t cmp_len = std::min(m, tlen);
        
        // Check if this suffix starts with pattern
        int c = std::char_traits<char>::compare(T + sa_pos, pattern.data(), cmp_len);
        if (c != 0 || cmp_len < m) {
            // Either doesn't match, or suffix is shorter than pattern
            break;
        }
        ++last_it;
    }

    for (auto it = first_it; it != last_it; ++it) {
        result.push_back(*it); // *it is a text position
    }

    return result;
}

// --- Return (text positions, SA indices) for all matches ---
std::pair<std::vector<uint64_t>, std::vector<uint64_t>>
SuffixArray::search_val_and_pos(const std::string& pattern) const {
    std::vector<uint64_t> values;
    std::vector<uint64_t> positions;
    if (pattern.empty() || suffixArray.empty()) return {values, positions};

    const char* T = text.data();
    const size_t N = text.size();
    const size_t m = pattern.size();

    // Find the first suffix that is >= pattern
    auto first_it = std::lower_bound(
        suffixArray.begin(), suffixArray.end(), pattern,
        [&](uint64_t sa_pos, const std::string& p) {
            size_t tlen = N - sa_pos;
            size_t cmp_len = std::min(m, tlen);
            int c = std::char_traits<char>::compare(T + sa_pos, p.data(), cmp_len);
            if (c != 0) return c < 0;
            return tlen < m;
        });

    // Find all suffixes that start with pattern
    auto last_it = first_it;
    while (last_it != suffixArray.end()) {
        uint64_t sa_pos = *last_it;
        size_t tlen = N - sa_pos;
        size_t cmp_len = std::min(m, tlen);
        
        int c = std::char_traits<char>::compare(T + sa_pos, pattern.data(), cmp_len);
        if (c != 0 || cmp_len < m) {
            break;
        }
        ++last_it;
    }

    for (auto it = first_it; it != last_it; ++it) {
        values.push_back(*it); // text position
        positions.push_back(static_cast<uint64_t>(
            static_cast<size_t>(it - suffixArray.begin()))); // SA index
    }
    
    return {values, positions};
}

void SuffixArray::save(const std::string& filename) const {
    std::ofstream out(filename, std::ios::binary);
    if (!out) {
        throw std::runtime_error("Cannot open file for writing: " + filename);
    }
    
    // 1. Save text
    size_t text_size = text.size();
    out.write(reinterpret_cast<const char*>(&text_size), sizeof(text_size));
    out.write(text.data(), text_size);
    
    // 2. Save suffix array
    size_t sa_size = suffixArray.size();
    out.write(reinterpret_cast<const char*>(&sa_size), sizeof(sa_size));
    out.write(reinterpret_cast<const char*>(suffixArray.data()), 
              sa_size * sizeof(uint64_t));
    
    // 3. Save LCP array
    size_t lcp_size = lcpArray.size();
    out.write(reinterpret_cast<const char*>(&lcp_size), sizeof(lcp_size));
    out.write(reinterpret_cast<const char*>(lcpArray.data()), 
              lcp_size * sizeof(int));
    
    out.close();
    
    if (!out.good()) {
        throw std::runtime_error("Error writing to file: " + filename);
    }
}

SuffixArray SuffixArray::load(const std::string& filename) {
    std::ifstream in(filename, std::ios::binary);
    if (!in) {
        throw std::runtime_error("Cannot open file for reading: " + filename);
    }
    
    SuffixArray result;
    
    // 1. Load text
    size_t text_size;
    in.read(reinterpret_cast<char*>(&text_size), sizeof(text_size));
    result.text.resize(text_size);
    in.read(&result.text[0], text_size);
    
    // 2. Load suffix array
    size_t sa_size;
    in.read(reinterpret_cast<char*>(&sa_size), sizeof(sa_size));
    result.suffixArray.resize(sa_size);
    in.read(reinterpret_cast<char*>(result.suffixArray.data()), 
            sa_size * sizeof(uint64_t));
    
    // 3. Load LCP array
    size_t lcp_size;
    in.read(reinterpret_cast<char*>(&lcp_size), sizeof(lcp_size));
    result.lcpArray.resize(lcp_size);
    in.read(reinterpret_cast<char*>(result.lcpArray.data()), 
            lcp_size * sizeof(int));
    
    in.close();
    
    if (!in.good()) {
        throw std::runtime_error("Error reading from file: " + filename);
    }
    
    return result;
}