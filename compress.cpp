#include "compress.hpp"
#include "suffix_array.hpp"
#include <iostream>
#include <algorithm>
#include <iterator>  


void Compress::initHash(const std::string& kmer)
{
    this ->hashMap[kmer] = hashValue();   
}

void Compress::initHashMap(SuffixArray& sa, const int k)
{
    const std::vector<int>& suffixArray = sa.getSuffixArray();
    std::vector<int>& lcpArray = sa.getLCPArray();
    const std::string& text = sa.getText();
    unsigned long i = 0;

    while(i < suffixArray.size())
    {
        int pos = suffixArray[i];
        std::string kmer_old = text.substr(pos, k);
        if (static_cast<int>(kmer_old.size()) < k || kmer_old.find('$') != std::string::npos)
        {
            std::cout << "Kmer is shorter than k: " << kmer_old << "\n";
            i++;
            continue;
        }

        

        for (unsigned long j = 0; j < kmer_old.size(); j++)
        {
            std::cout << kmer_old[j];
        }
        std::cout << "\n";

        unsigned long steps = i;
        // std::cout << "i: " << i ;
        Compress::initHash(kmer_old);
        while (lcpArray[steps] >= k) steps++;
        // std::cout << " step: " << step << "\n";
        
        
        i = steps + 1;
        

        
    }
    std::cout << "Map size: " << hashMap.size() << "\n";
}

void Compress::printMap(){
    for (const auto& [key, value] : this -> hashMap) {
        std::cout << "Key: \"" << key << "\""
                  << ", kSAindex: " << value.kSAindex
                  << ", occurences: " << value.occurences
                  << ", shift: " << value.shift
                  << ", trace: " << value.traceback 
                  << std::endl;
    }
}

unsigned long Compress::findLCPmaxmin(SuffixArray& sa,const int k)
{
    std::vector<int>& lcpArray = sa.getLCPArray();
    auto maxIt = std::max_element(lcpArray.begin(), lcpArray.end());
    int index = std::distance(lcpArray.begin(), maxIt);

    unsigned long upperBound = index;
    unsigned long lowerBound = index;

    while (lcpArray[upperBound + 1] >= k && upperBound < lcpArray.size() - 2) {
        upperBound++;
    }

    while (lcpArray[lowerBound - 1] >= k && lowerBound > 1) {
        lowerBound--;
    }

    auto minIt = std::min_element(lcpArray.begin() + lowerBound, lcpArray.begin() + upperBound);
    unsigned long minIndex = std::distance(lcpArray.begin(), minIt);
    std::fill(lcpArray.begin() + lowerBound, lcpArray.begin() + upperBound + 2, -1);
    sa.printLCPArray();

    return minIndex;
    
}

void Compress::compression (SuffixArray& sa, const int k)
{
    Compress cSA;
    cSA.initHashMap(sa, k);
    sa.printSuffixArray();
    const std::vector<int>& SuffixArray = sa.getSuffixArray();
    const std::string& text = sa.getText();

    unsigned long startIndex = Compress::findLCPmaxmin(sa, k);
    std::cout << "Minimum LCP index for k = " << k << ": " << startIndex << "\n";

    std::string pattern = text.substr(SuffixArray[startIndex]);
    std::cout << "Pattern: " << pattern << "\n";

    cSA.readPattern(sa, cSA, pattern, k);

    cSA.printMap();

    







}

void Compress::readPattern (SuffixArray& sa, Compress& cSA, std::string& pattern, int k)
{
    std::string kmer_old = pattern.substr(0, k);
    std::cout << "Kmer_old: " << kmer_old << "\n";
    std::pair<std::vector<int>,std::vector<int>> val_and_index = sa.search_index_and_value(kmer_old);
    std::vector<int> tempResult = val_and_index.first;
    std::vector<int> tempIndex = val_and_index.second;
    cSA.compressedSA = tempResult;
    // std::cout << "Result size: " << tempResult.size() << "\n";
    cSA.hashMap[kmer_old].kSAindex = 0;
    cSA.hashMap[kmer_old].occurences = tempResult.size();
    cSA.hashMap[kmer_old].shift = -1; // no shift, because it is not dependent on another pattern

    int shift = 1;

    unsigned old_res_size = cSA.compressedSA.size();
    std::vector<int> lcpArray = sa.getLCPArray();
    
    for (unsigned i = 1; i < pattern.size() - 2;i++)
    {
        std::string kmer_new = pattern.substr(i,k);
        std::cout << "Kmer: " << kmer_new << "\n";
        val_and_index = sa.search_index_and_value(kmer_new);
        tempResult = val_and_index.first;
        tempIndex = val_and_index.second;
        unsigned i_res = 0;
        unsigned i_temp = 0;

        while (i_res < old_res_size && i_temp < tempResult.size())
        {   
            //   std::cout << "Flag: " << i << "\n";
            //     std::cout << "cSA.compressedSA[i_res] + shift : " << cSA.compressedSA[i_res] + shift << "\n";
            //     std::cout << "tempResult[i_temp]: " << tempResult[i_temp] << "\n";
            //     std::cout << "shift: " << shift << "\n";
            if (cSA.compressedSA[i_res] + shift == tempResult[i_temp] )
            {   
                lcpArray[i_temp] = -1;
                i_res++;
                i_temp++;
                
            }
            else
            {   
                cSA.compressedSA.push_back(tempResult[i_temp]);
                cSA.hashMap[kmer_new].kSAindex = cSA.compressedSA.size() - 1;
                cSA.hashMap[kmer_new].occurences += 1;
                lcpArray[i_temp] = -1;
                i_temp++;
            }
            sa.printLCPArray();
           
        }
        std::cout << "i_res: " << i_res << "\n";
        std::cout << "old_res_size: " << old_res_size << "\n";
        if(i_res == old_res_size )
        {
            if(cSA.hashMap[kmer_old].traceback == nullptr)
            {   
                cSA.hashMap[kmer_new].traceback = &cSA.hashMap[kmer_old];
            }
            else
            {
                cSA.hashMap[kmer_new].traceback = cSA.hashMap[kmer_old].traceback;
            }
        }

        while (i_temp < tempResult.size())
        {
            cSA.compressedSA.push_back(tempResult[i_temp]);
            cSA.hashMap[kmer_new].kSAindex = cSA.compressedSA.size() - 1;
            cSA.hashMap[kmer_new].occurences += 1;
            i_temp++;
        }
        cSA.hashMap[kmer_new].shift = shift;
        kmer_old = kmer_new;
        shift++;
    }
}

int main()
{
    SuffixArray sa("ATAACCGA$ATGACCGA$ATAACCGA$CTAACCGA$ATAACCGA");
    Compress c = Compress();
    c.initHashMap(sa, 3);

    c.compression(sa, 3);
    return 0;


}