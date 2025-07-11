#include <iostream>
#include "suffix_array.hpp"

int main() {
    std::string text;
    std::cout << "Geben Sie den Text ein: ";
    std::getline(std::cin, text);
    
    if (text.empty()) {
        std::cout << "Leerer Text eingegeben!\n";
        return 1;
    }
    
    // Erstelle Suffix Array
    SuffixArray sa(text);
    
    // Zeige Suffix Array
    sa.printSuffixArray();
    
    // Zeige LCP Array
    sa.printLCPArray();
    
    // Beispiel für Pattern-Suche
    std::string pattern;
    std::cout << "\nGeben Sie ein Pattern zum Suchen ein (oder Enter zum Beenden): ";
    std::getline(std::cin, pattern);
    
    if (!pattern.empty()) {
        std::vector<int> matches = sa.search(pattern);
        
        if (matches.empty()) {
            std::cout << "Pattern '" << pattern << "' nicht gefunden.\n";
        } else {
            std::cout << "Pattern '" << pattern << "' gefunden an Positionen: ";
            for (int pos : matches) {
                std::cout << pos << " ";
            }
            std::cout << "\n";
        }
    }
    
    return 0;
}