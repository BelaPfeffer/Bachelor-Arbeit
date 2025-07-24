#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <map>
#include "suffix_array.hpp"

    
    // Hilfsfunktion für das Sortieren der Suffixe
    void SuffixArray::buildSuffixArray() {
        int n = text.length();
        suffixArray.resize(n);
        
        // Initialisiere das Suffix Array mit Indizes
        for (int i = 0; i < n; i++) {
            suffixArray[i] = i;
        }
        
        // Sortiere die Suffixe lexikographisch
        std::sort(suffixArray.begin(), suffixArray.end(), [this](int a, int b) {
            return text.substr(a) < text.substr(b);
        });
    }
    
    size_t SuffixArray::memoryUsageBytes() const {
        size_t size = sizeof(*this);

        size += text.capacity() * sizeof(char);
        size += suffixArray.capacity() * sizeof(int);
        size += lcpArray.capacity() * sizeof(int);

        return size;
    }

    // Effizientere Implementierung mit O(n log n) Komplexität
    void SuffixArray::buildSuffixArrayOptimized() {
        int n = text.length();
        suffixArray.resize(n);
        
        // Rank Array für die aktuelle Iteration
        std::vector<int> rank(n);
        std::vector<int> tempRank(n);
        
        // Initialisiere Ranks basierend auf einzelnen Zeichen
        for (int i = 0; i < n; i++) {
            rank[i] = text[i];
            suffixArray[i] = i;
        }
        
        // Verdopple die Länge der zu vergleichenden Substrings
        for (int k = 1; k < n; k *= 2) {
            // Sortiere basierend auf dem aktuellen Rank
            std::sort(suffixArray.begin(), suffixArray.end(), [&](int a, int b) {
                if (rank[a] != rank[b]) {
                    return rank[a] < rank[b];
                }
                int rankA = (a + k < n) ? rank[a + k] : -1;
                int rankB = (b + k < n) ? rank[b + k] : -1;
                return rankA < rankB;
            });
            
            // Berechne neue Ranks
            tempRank[suffixArray[0]] = 0;
            for (int i = 1; i < n; i++) {
                tempRank[suffixArray[i]] = tempRank[suffixArray[i-1]];
                
                int prev = suffixArray[i-1];
                int curr = suffixArray[i];
                
                if (rank[prev] != rank[curr] || 
                    (prev + k < n ? rank[prev + k] : -1) != (curr + k < n ? rank[curr + k] : -1)) {
                    tempRank[suffixArray[i]]++;
                }
            }
            
            rank = tempRank;
        }
    }
    
    // Berechne das LCP Array
    void SuffixArray::buildLCPArray() {
        int n = text.length();
        lcpArray.resize(n - 1);
        
        // Erstelle inverse Suffix Array für O(1) Zugriff
        std::vector<int> invSuffixArray(n);
        for (int i = 0; i < n; i++) {
            invSuffixArray[suffixArray[i]] = i;
        }
        
        int lcp = 0;
        for (int i = 0; i < n; i++) {
            if (invSuffixArray[i] == n - 1) {
                lcp = 0;
                continue;
            }
            
            int j = suffixArray[invSuffixArray[i] + 1];
            
            while (i + lcp < n && j + lcp < n && text[i + lcp] == text[j + lcp]) {
                lcp++;
            }
            
            lcpArray[invSuffixArray[i]] = lcp;
            
            if (lcp > 0) {
                lcp--;
            }
        }
    }
    
    // Konstruktor
    SuffixArray::SuffixArray(const std::string& input) : text(input) {
        if (!text.empty()) {
            buildSuffixArrayOptimized();
            buildLCPArray();
        }
    }
    
    // Getter für das Suffix Array
    const std::vector<int>& SuffixArray::getSuffixArray() const {
        return suffixArray;
    }
    const std::string& SuffixArray::getText() const {
        return text;
    }
    
    // Getter für das LCP Array
    const std::vector<int>& SuffixArray::getLCPArray() const {
        return lcpArray;
    }
    
    // Ausgabe des Suffix Arrays
    void SuffixArray::printSuffixArray() const {
        std::cout << "Suffix Array:\n";
        std::cout << "Index\tPosText\tLCP\tSuffix\n";
        for (unsigned long i = 0; i < suffixArray.size(); i++) {
            std::cout << i << "\t"<<suffixArray[i] << "\t" << lcpArray[i] << "\t" << text.substr(suffixArray[i]) << "\n";
        }
    }
    
    // Ausgabe des LCP Arrays
    void SuffixArray::printLCPArray() const {
        std::cout << "\nLCP Array:\n";
        std::cout << "Index\tLCP\tSuffix 1\t\tSuffix 2\n";
        for (unsigned long i = 0; i < lcpArray.size(); i++) {
            std::cout << i << "\t" << lcpArray[i] << "\t" 
                      << text.substr(suffixArray[i]) << "\t\t" 
                      << text.substr(suffixArray[i + 1]) << "\n";
        }
    }
    
    // Suche nach einem Pattern (binäre Suche)
    std::vector<int> SuffixArray::search(const std::string& pattern) const {
        std::vector<int> result;
        int n = suffixArray.size();
        int patternLen = pattern.length();
        
        // Binäre Suche für den ersten Treffer
        int left = 0, right = n - 1;
        int first = -1;
        
        while (left <= right) {
            int mid = (left + right) / 2;
            std::string suffix = text.substr(suffixArray[mid], patternLen);
            
            if (suffix >= pattern) {
                if (suffix == pattern) {
                    first = mid;
                }
                right = mid - 1;
            } else {
                left = mid + 1;
            }
        }
        
        if (first == -1) {
            return result;
        }
        
        // Binäre Suche für den letzten Treffer
        left = 0;
        right = n - 1;
        int last = -1;
        
        while (left <= right) {
            int mid = (left + right) / 2;
            std::string suffix = text.substr(suffixArray[mid], patternLen);
            
            if (suffix <= pattern) {
                if (suffix == pattern) {
                    last = mid;
                }
                left = mid + 1;
            } else {
                right = mid - 1;
            }
        }
        
        // Sammle alle Treffer
        for (int i = first; i <= last; i++) {
            result.push_back(suffixArray[i]);
        }
        
        return result;
    }


std::pair<std::vector<int>, std::vector<int>> SuffixArray::search_val_and_pos(const std::string& pattern) const {
        std::vector<int> values;
        std::vector<int> positions;
        int n = suffixArray.size();
        int patternLen = pattern.length();
        
        // Binäre Suche für den ersten Treffer
        int left = 0, right = n - 1;
        int first = -1;
        
        while (left <= right) {
            int mid = (left + right) / 2;
            std::string suffix = text.substr(suffixArray[mid], patternLen);
            
            if (suffix >= pattern) {
                if (suffix == pattern) {
                    first = mid;
                }
                right = mid - 1;
            } else {
                left = mid + 1;
            }
        }
        
        if (first == -1) {
            return {values, positions};
        }
        
        // Binäre Suche für den letzten Treffer
        left = 0;
        right = n - 1;
        int last = -1;
        
        while (left <= right) {
            int mid = (left + right) / 2;
            std::string suffix = text.substr(suffixArray[mid], patternLen);
            
            if (suffix <= pattern) {
                if (suffix == pattern) {
                    last = mid;
                }
                left = mid + 1;
            } else {
                right = mid - 1;
            }
        }
        
        // Sammle alle Treffer
        for (int i = first; i <= last; i++) {
            values.push_back(suffixArray[i]);
            positions.push_back(i);
        }

        return {values, positions};
    }