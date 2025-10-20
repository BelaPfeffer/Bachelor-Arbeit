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

// --- Binary search over SA without creating substrings ---
std::vector<uint64_t> SuffixArray::search(const std::string& pattern) const {
    std::vector<uint64_t> result;
    if (pattern.empty() || suffixArray.empty()) return result;

    const char* T = text.data();
    const size_t N = text.size();
    const size_t m = pattern.size();

    auto suffix_less_than_pattern = [&](uint64_t sa_pos, const std::string& p) {
        const size_t tlen  = N - static_cast<size_t>(sa_pos);
        const size_t limit = (m < tlen) ? m : tlen;
        const int c = std::char_traits<char>::compare(T + sa_pos, p.data(), limit);
        if (c != 0) return c < 0;
        // equal up to min; shorter string is lexicographically smaller
        return tlen < m;
    };

    auto pattern_less_than_suffix = [&](const std::string& p, uint64_t sa_pos) {
        const size_t tlen  = N - static_cast<size_t>(sa_pos);
        const size_t limit = (m < tlen) ? m : tlen;
        const int c = std::char_traits<char>::compare(p.data(), T + sa_pos, limit);
        if (c != 0) return c < 0;
        return m < tlen;
    };

    auto first_it = std::lower_bound(
        suffixArray.begin(), suffixArray.end(), pattern, suffix_less_than_pattern);

    auto last_it  = std::upper_bound(
        first_it, suffixArray.end(), pattern, pattern_less_than_suffix);

    for (auto it = first_it; it != last_it; ++it)
        result.push_back(*it); // *it is a text position

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

    auto suffix_less_than_pattern = [&](uint64_t sa_pos, const std::string& p) {
        const size_t tlen  = N - static_cast<size_t>(sa_pos);
        const size_t limit = (m < tlen) ? m : tlen;
        const int c = std::char_traits<char>::compare(T + sa_pos, p.data(), limit);
        if (c != 0) return c < 0;
        return tlen < m;
    };

    auto pattern_less_than_suffix = [&](const std::string& p, uint64_t sa_pos) {
        const size_t tlen  = N - static_cast<size_t>(sa_pos);
        const size_t limit = (m < tlen) ? m : tlen;
        const int c = std::char_traits<char>::compare(p.data(), T + sa_pos, limit);
        if (c != 0) return c < 0;
        return m < tlen;
    };

    auto first_it = std::lower_bound(
        suffixArray.begin(), suffixArray.end(), pattern, suffix_less_than_pattern);

    auto last_it  = std::upper_bound(
        first_it, suffixArray.end(), pattern, pattern_less_than_suffix);

    for (auto it = first_it; it != last_it; ++it) {
        // Optional equality guard
        if (text.compare(static_cast<size_t>(*it), m, pattern) == 0) {
            values.push_back(*it); // text position
            positions.push_back(static_cast<uint64_t>(
                static_cast<size_t>(it - suffixArray.begin()))); // SA index
        }
    }
    return {values, positions};
}
