#include "compress.hpp"
#include "suffix_array.hpp"
#include <iostream>
#include <algorithm>
#include <iterator>  

void CompressedSA::printCompressedSA()
{
    std::cout << "Compressed Suffix Array:\n";
    for (unsigned long i = 0; i < compressedSA.size(); i++) {
        std::cout << i << "\t" << compressedSA[i] << "\t" << this -> text.substr(compressedSA[i]) << "\n";
    }
    std::cout << "\n";
}

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

        if (pos + k > text.size() )
        {
            // std::cout << "Kmer is shorter than k: " << kmer_old << "\n";
            i++;
            continue;
        }

         std::string kmer = text.substr(pos, k);

        if(kmer.find('$') != std::string::npos)
        {
            i++;
            continue;
        }

    
        for (unsigned long j = 0; j < kmer.size(); j++)
        {
            std::cout << kmer_old[j];
        }
        std::cout << "\n";

        unsigned long steps = i;
        // std::cout << "i: " << i ;
        this -> initHash(kmer_old);

        //jumps to next Bucket
        while (lcpArray[steps] >= k) steps++;
        
        
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
                  << ", processed: " << value.processed
                  << std::endl;
    }
}

std::pair<unsigned long, int> CompressedSA::findLCPmaxmin(const int k)
{
    std::vector<int>& lcpArray = this -> lcpArray;
    auto maxIt = std::max_element(lcpArray.begin(), lcpArray.end());
    int index = std::distance(lcpArray.begin(), maxIt);

    unsigned long upperBound = index;
    unsigned long lowerBound = index;

    while (lcpArray[upperBound + 1] >= k && upperBound <= lcpArray.size() - k) {
        upperBound++;
    }

    while (lcpArray[lowerBound - 1] >= k && lowerBound > 1) {
        lowerBound--;
    }


    auto minIt = std::min_element(lcpArray.begin() + lowerBound, lcpArray.begin() + upperBound);
    unsigned long minIndex = std::distance(lcpArray.begin(), minIt);
    std::pair<unsigned long, int>minMax = {minIndex, lcpArray[minIndex]};
    std::fill(lcpArray.begin() + lowerBound, lcpArray.begin() + upperBound + k-1, -1);
    // this -> printLCPArray();

    return minMax;
    
}

void CompressedSA::compression (const int k)
{
    this -> initHashMap(k);
    // this -> printSuffixArray();
    const std::vector<int>& SuffixArray = this -> suffixArray;
    const std::string& text = this -> text;
    std::pair<unsigned long, int> start = findLCPmaxmin (k);
    unsigned long startIndex = start.first;
    int LCPVal = start.second;
    this -> printLCPArray();

    std::cout << "LCPVal: " << LCPVal << "\n";

    while (LCPVal > -1)
    {
       
        std::cout << "Minimum LCP index for k = " << k << ": " << startIndex << "\n";
        std::string pattern = text.substr(SuffixArray[startIndex]);
        if (pattern.find('$') != std::string::npos || static_cast<int>(pattern.size()) < k)
        {
            
            std::cout << "Pattern contains end character '$', skipping...\n";
            this -> lcpArray[startIndex] = -1; // Mark as processed
            start = findLCPmaxmin(k);
            startIndex = start.first;
            LCPVal = start.second;
            std::cout << "Next startIndex: " << startIndex << "\n";
            continue;
        }
        std::cout << "Pattern: " << pattern << "\n";

        this -> readPattern(pattern, k);

        this -> printMap();
        
        this -> printLCPArray();

        start = findLCPmaxmin(k);
        startIndex = start.first;
        LCPVal = start.second;
    }

    std::cout << "Compression finished.\n";
    this -> printMap();
    this -> printLCPArray();
    this -> printCompressedSA();
}

void CompressedSA::readPattern (std::string& pattern, int k)
{
    std::string kmer_old = pattern.substr(0, k);
    std::cout << "Kmer_old: " << kmer_old << "\n";
    std::pair<std::vector<int>, std::vector<int>> tempRes = this -> search_val_and_pos(kmer_old);
    this -> compressedSA.insert(compressedSA.end(), tempRes.first.begin(), tempRes.first.end());
    std::cout << "CompressedSA: "   << "\n";
    this -> printCompressedSA();
    // std::cout << "Result size: " << tempResult.size() << "\n";
    this -> hashMap[kmer_old].kSAindex = 0;
    this -> hashMap[kmer_old].occurences = tempRes.first.size();
    this -> hashMap[kmer_old].shift = -1; // no shift, because it is not dependent on another pattern
    this -> hashMap[kmer_old].processed = true;

    int shift = 1;

    unsigned old_res_size = this -> compressedSA.size();
    
    for (unsigned i = 1; i <= pattern.size() - k;i++)
    {
        std::string kmer_new = pattern.substr(i,k);
        std::cout << "Kmer: " << kmer_new << "\n";

        if (this -> hashMap[kmer_new].processed)
        {
            std::cout << "Kmer already processed, skipping...\n";
            break;
        }

        tempRes = this -> search_val_and_pos(kmer_new);
        unsigned i_res = 0;
        unsigned i_temp = 0;

        while (i_res < old_res_size && i_temp < tempRes.first.size())
        {   
            //   std::cout << "Flag: " << i << "\n";
                std::cout << "cSA.compressedSA[i_res] + shift : " << this -> compressedSA[i_res] + shift << "\n";
                std::cout << "tempRes.first[i_temp]: " << tempRes.first[i_temp] << "\n";
                std::cout << "shift: " << shift << "\n";
            if (this -> compressedSA[i_res] + shift == tempRes.first[i_temp] )
            {   
                std::cout << "Chaining to compressedSA: " << this -> compressedSA[i_res] << "\n";
                this -> lcpArray[tempRes.second[i_temp]] = -1;
                i_res++;
                i_temp++;
            }
            else
            {   
                std::cout << "Adding to compressedSA: " << tempRes.first[i_temp] << "\n";
                this -> compressedSA.push_back(tempRes.first[i_temp]);
                this -> hashMap[kmer_new].kSAindex = this -> compressedSA.size() - 1;
                this -> hashMap[kmer_new].occurences += 1;
                this -> lcpArray[tempRes.second[i_temp]] = -1;
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

        while (i_temp < tempRes.first.size())
        {
            this -> compressedSA.push_back(tempRes.first[i_temp]);
            this -> hashMap[kmer_new].kSAindex = this -> compressedSA.size() - 1;
            this -> hashMap[kmer_new].occurences += 1;
            i_temp++;
        }
        this -> hashMap[kmer_new].shift = shift;
        this -> hashMap[kmer_new].processed = true;
        kmer_old = kmer_new;
        shift++;
    }
}

int main()
{
    CompressedSA cSA("ATAACCGA$ATGACCGA$ATAACCGA$CTAACCGA$ATAACCGA");
    // CompressedSA cSA("ATAACCGA$ATGACCGA");
    // CompressedSA c = CompressedSA();
    // cSA.initHashMap(3);
   

    cSA.compression(3);
    
    return 0;


}