#include "compress.hpp"
#include "suffix_array.hpp"
#include <iostream>
#include <algorithm>
#include <iterator>  


void CompressedSA::initHash(const std::string& kmer)
{
    this ->hashMap[kmer] = hashValue();   
}

void CompressedSA::initHashMap(const int k)
{
    const std::vector<int>& suffixArray = this->suffixArray;
    const std::vector<int>& lcpArray = this->lcpArray;
    const std::string& text = this->text;
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
        this -> initHash(kmer_old);
        while (lcpArray[steps] >= k) steps++;
        // std::cout << " step: " << step << "\n";
        
        
        i = steps + 1;
        

        
    }
    std::cout << "Map size: " << hashMap.size() << "\n";
}

void CompressedSA::printMap(){
    for (const auto& [key, value] : this -> hashMap) {
        std::cout << "Key: \"" << key << "\""
                  << ", kSAindex: " << value.kSAindex
                  << ", occurences: " << value.occurences
                  << ", shift: " << value.shift
                  << ", trace: " << value.traceback 
                  << std::endl;
    }
}

unsigned long CompressedSA::findLCPmaxmin(const int k)
{
    const std::vector<int>& lcpArray = this -> lcpArray;
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

    return minIndex;
    
}

void CompressedSA::compression (const int k)
{
    this -> initHashMap(k);
    this -> printSuffixArray();
    const std::vector<int>& SuffixArray = this -> suffixArray;
    const std::string& text = this -> text;

    unsigned long startIndex = findLCPmaxmin (k);
    std::cout << "Minimum LCP index for k = " << k << ": " << startIndex << "\n";

    std::string pattern = text.substr(SuffixArray[startIndex]);
    std::cout << "Pattern: " << pattern << "\n";

    this -> readPattern(pattern, k);

    this -> printMap();

    







}

void CompressedSA::readPattern (std::string& pattern, int k)
{
    std::string kmer_old = pattern.substr(0, k);
    std::cout << "Kmer_old: " << kmer_old << "\n";
    std::vector<int> tempResult = this -> search(kmer_old);
    this -> compressedSA = tempResult;
    // std::cout << "Result size: " << tempResult.size() << "\n";
    this -> hashMap[kmer_old].kSAindex = 0;
    this -> hashMap[kmer_old].occurences = tempResult.size();
    this -> hashMap[kmer_old].shift = -1; // no shift, because it is not dependent on another pattern

    int shift = 1;

    unsigned old_res_size = this -> compressedSA.size();
    
    for (unsigned i = 1; i <= pattern.size() - k;i++)
    {
        std::string kmer_new = pattern.substr(i,k);
        std::cout << "Kmer: " << kmer_new << "\n";
        tempResult = this -> search(kmer_new);
        unsigned i_res = 0;
        unsigned i_temp = 0;

        while (i_res < old_res_size && i_temp < tempResult.size())
        {   
              std::cout << "Flag: " << i << "\n";
                std::cout << "cSA.compressedSA[i_res] + shift : " << this -> compressedSA[i_res] + shift << "\n";
                std::cout << "tempResult[i_temp]: " << tempResult[i_temp] << "\n";
                std::cout << "shift: " << shift << "\n";
            if (this ->compressedSA[i_res] + shift == tempResult[i_temp] )
            {   
              
                i_res++;
                i_temp++;
            }
            else
            {   
                this -> compressedSA.push_back(tempResult[i_temp]);
                this -> hashMap[kmer_new].kSAindex = this -> compressedSA.size() - 1;
                this -> hashMap[kmer_new].occurences += 1;
                i_temp++;
            }
           
        }
        
        if(i_res == old_res_size )
        {
            if(this -> hashMap[kmer_old].traceback == nullptr)
            {   
                this -> hashMap[kmer_new].traceback = &this -> hashMap[kmer_old];
            }
            else
            {
                this -> hashMap[kmer_new].traceback = this ->hashMap[kmer_old].traceback;
            }
        }

        while (i_temp < tempResult.size())
        {
            this -> compressedSA.push_back(tempResult[i_temp]);
            this -> hashMap[kmer_new].kSAindex = this -> compressedSA.size() - 1;
            this -> hashMap[kmer_new].occurences += 1;
            i_temp++;
        }
        this -> hashMap[kmer_new].shift = shift;
        kmer_old = kmer_new;
        shift++;
    }
}

int main()
{
    CompressedSA cSA("ATAACCGA$ATGACCGA$ATAACCGA$CTAACCGA$ATAACCGA");
    // CompressedSA c = CompressedSA();
    cSA.initHashMap(3);
    cSA.printMap();

    cSA.compression(3);
    return 0;


}