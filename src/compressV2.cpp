#include "compressV2.hpp"
#include "suffix_array.hpp"
#include <iostream>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <iterator>  

#include <sdsl/suffix_array_algorithm.hpp>

void CompressedSA::setBits(lcp_interval& iv, bool bitVal)
{
    std::fill(computeSuffix.begin() + iv.left, computeSuffix.begin() + iv.right + 1, bitVal);
    sdsl::util::init_support(rankSupport, &computeSuffix); // Reinitialize rank support after modifying computeSuffix
}

void CompressedSA::printCompressedSA()
{
    std::cout << "Compressed Suffix Array:\n";
    for (unsigned long i = 0; i < compressedSA.size(); i++) {
        std::cout << i << "\t" << compressedSA[i] << "\t" << this -> text.substr(compressedSA[i]) << "\n";
    }
    std::cout << "\n";
}

 size_t CompressedSA::memoryUsageBytes() const {
    size_t size = sizeof(*this); // Speicher für CompressedSA selbst

    // Dynamisch belegter Speicher im Vektor
    size += compressedSA.capacity() * sizeof(int);

    // Dynamisch belegter Speicher in der Hashmap
    size += hashMap.size() * (sizeof(uint64_t) + sizeof(hashValue));

    // Optional: Wenn load factor hoch ist, ist capacity > size:
    size += (hashMap.bucket_count() - hashMap.size()) * (sizeof(void*) + sizeof(size_t)); // ungefähr

    // Speicher aus der Basisklasse (SuffixArray)
    size += SuffixArray::memoryUsageBytes(); 

    return size;
}


void CompressedSA::printMap(uint64_t k)
{
    for (const auto& [key, value] : this -> hashMap) {

        // std::string csaIndexes = "[";
        // for (const auto& index : value.cSAindex) {
        //     csaIndexes += std::to_string(index) + ",";
        // }

        // csaIndexes += "]";

        // if (value.processed == false) {
        //     continue; // Skip uninitialized entries
        // }
        std::cout << "Key: \"" << key << "\""
                  << ", Decoded_Key: " << decode_dna5(key,k) // Assuming k=2 for decoding
                  << ", cSAindex: " << value.cSAindex
                  << ", occurences: " << value.occurences
                  << ", lcp_interval_index: " << value.lcp_interval_index
                  << ", shift: " << value.shift
                  << ", refOcc: " << value.refOccurrences
                  << ", trace: " << value.traceback_key 
                  << ", processed: " << value.processed
                  << std::endl;
    }
}

void CompressedSA::initHash(const uint64_t& kmer, unsigned long lcp_interval_index)
{
    this ->hashMap[kmer] = hashValue();
    this ->hashMap[kmer].lcp_interval_index = lcp_interval_index;

}

// void CompressedSA::initHashMap(const unsigned k)
// {
//     unsigned long i = 0;

//     while(i < suffixArray.size())
//     {
//         unsigned long pos = suffixArray[i];

//         if (pos + k > text.size() )
//         {
//             // std::cout << "Kmer is shorter than k: " << kmer_old << "\n";
//             i++;
//             continue;
//         }

//         std::string kmer = text.substr(pos, k);
//         // std::cout << "Kmer: " << kmer << "\n";
//         // std::cout << "pos: " << pos << "\n";

//         if(kmer.find('$') != std::string::npos)
//         {
//             i++;
//             continue;
//         }

//         uint64_t encodedKmer = encode_dna5(kmer);
//         this -> initHash(encodedKmer, lcpIntervals.size() - 1);
//         for (unsigned long j = 0; j < kmer.size(); j++)
//         {
//             std::cout << kmer[j];
//         }
//         std::cout << "\n";

//         unsigned long steps = i + 1;
//         if(steps == suffixArray.size())
//         {
//             break;
//         }
//         // std::cout << "steps: " << steps << "\n";
        

        
//         while (static_cast<unsigned>(lcpArray[steps]) >= k) {steps++; };
//         // std::cout << "lcp loopend:" << std::endl;
   


//         i = steps;

//     }
//     std::cout << "Map size: " << hashMap.size() << "\n";
// }

uint64_t CompressedSA::encode_dna5(const std::string& kmer)
{
    uint64_t encoded = 0;
    for (char c : kmer) {
        encoded <<= 3; // Shift left by  bits
        switch (c) {
            case 'A': encoded |= 0; break; // 000
            case 'C': encoded |= 1; break; // 001
            case 'G': encoded |= 2; break; // 010
            case 'T': encoded |= 3; break; // 011
            case '$': encoded |= 4; break; // 100 (special end character)
            default: throw std::invalid_argument("Invalid character in DNA sequence");
        }
    }
    return encoded;
}

std::string CompressedSA::decode_dna5(uint64_t encoded, unsigned k)
{
    std::string decoded;
    for (int i = 0; i < k; ++i) {
        int base = encoded & 0x07; // Get the last 3 bits
        switch (base) {
            case 0: decoded += 'A'; break;
            case 1: decoded += 'C'; break;
            case 2: decoded += 'G'; break;
            case 3: decoded += 'T'; break;
            case 4: decoded += '$'; break; // Special end character
            default: throw std::invalid_argument("Invalid base in encoded DNA sequence");
        }
        encoded >>= 3; // Shift right by 3 bits
    }
    std::reverse(decoded.begin(), decoded.end()); // Reverse to get the correct order
    return decoded;
}

hashValue CompressedSA::getHashMapValue(const std::string& kmer)
{
    uint64_t encodedKmer = encode_dna5(kmer);
    return this->hashMap[encodedKmer];
}

void CompressedSA::setReferenceValue(const uint64_t encodedKmer, const int shift, unsigned long refOcurrences, const uint64_t traceback_key, const bool processed)
{
    this -> hashMap[encodedKmer].setReferenceValue(shift, refOcurrences, traceback_key, processed);
}


void CompressedSA::setValue(const uint64_t encodedKmer, const unsigned long cSAindex, const unsigned long occurences) // Index im komprimierten Suffix Array
{
   this -> hashMap[encodedKmer].setValue(cSAindex, occurences);
}

void CompressedSA::printSuffixArray() const 
{
        std::cout << "Suffix Array:\n";
        std::cout << "Index\tPosText\tLCP\tbitvektor\tSuffix\n";
        for (unsigned long i = 0; i < suffixArray.size(); i++) {
            std::cout << i << "\t"<<suffixArray[i] << "\t" << lcpArray[i] << "\t" << computeSuffix[i] << "\t" << text.substr(suffixArray[i]) << "\n";
        }
}

lcp_interval CompressedSA::get_lcp_interval(unsigned i, unsigned k)
{
    const sdsl::lcp_bitcompressed<>& lcp = lcpArray;
    unsigned long n = lcp.size() + 1; // Because LCP array has size n - 1
    unsigned long left = i;
    unsigned long right = i;
    unsigned long min_index = i;
    
    // Expand to the left
    while (left > 0 && lcp[left] >= k) {
        --left;
        if (lcp[left] < lcp[min_index]) {
            min_index = left;
        }
    }
    
    // Expand to the right
    while (right < n - 2 && lcp[right + 1] >= k) {
        ++right;
        if (lcp[right] < lcp[min_index]) {
            min_index = right;
        }
    }
    
    return lcp_interval(left, right, min_index);
}

//wird true wenn wert valid
bool CompressedSA::filterSA(const int k, unsigned i)
{
    unsigned long suffix = suffixArray[i];
    if (suffix + k > text.size()) {
        return false; // Suffix is too short to contain a k-mer of length k
    }
    else if( text.substr(suffix, k).find('$') != std::string::npos) {
        return false; // Suffix contains the end character '$'
    }
    return true;
}

void CompressedSA::initComputeSuffix(unsigned k)
{
    computeSuffix = bit_vector(suffixArray.size());
    for (unsigned long i = 0; i < suffixArray.size(); i++) {
        if (filterSA(k, i)) {
            computeSuffix[i] = 1; // Mark this suffix as valid
        } else {
            computeSuffix[i] = 0; // Mark this suffix as invalid
        }
    }
}

void CompressedSA::compression(const unsigned k, lcp_interval& interval)
{
    std::cout << "rankCalc: " << rankSupport(interval.right + 1) - rankSupport(interval.left) << "\n";
    // Check 
    if (rankSupport(interval.right + 1) - rankSupport(interval.left) == 0)
    {
        std::cout << "Interval already computed or no valid kmer in interval.\n";
        return;
    }

    //copies the suffixArray data in the compressedSA
    std::copy(suffixArray.begin() + interval.left, suffixArray.begin() + interval.right + 1, std::back_inserter(compressedSA));
    
    unsigned pat_pos_index = suffixArray[interval.min_index]; // Index des Musters im text
    unsigned long occurences = interval.right - interval.left + 1;  // Anzahl der Vorkommen des Musters im SuffixArray
    unsigned long current_csa_index = compressedSA.size() - occurences;


    uint64_t kmer_old = encode_dna5(text.substr(pat_pos_index, k));


    setValue(kmer_old, current_csa_index, occurences); // das könnte probleme geben wenn compressed SA noch leer ist
    lcpIntervals[hashMap[kmer_old].lcp_interval_index] = std::nullopt; 
    //wenn ich diese Zeilen Rausnehme könnte man Ringschluss bilden sodass die anzahl der Gespeicherten Suffixe noch kleiner wird. 
    //BSP.: Ich habe ACCGA als erstes Kompressionspattern benutzt später finde ich dann ATACCGA -> ACCGA kann in Abhängigkeit von 
    //ATACCGA dargestellt werden -> alle Pattern die in Abhängigkeit von ACCGA dargestellt werden können auch mit ATACCGA dargestellt werden.
    //Dazu müsste man nur ACCGA aus SA löschen und in Abhängigkeit von ATACCGA's Position im Text darstellen.

    setBits(interval, 0); // Markiere alle Suffixe im Intervall als computed

    // std::cout << "Suffix: " << text.substr(pat_pos_index) << "\n";

    // this -> printCompressedSA();
    // this -> printMap(k);

    // ab hier nochmal gucken

    int shift = 1;
    int mask = 0b100100100;

    unsigned text_index = pat_pos_index + shift;

    for (unsigned long i = pat_pos_index + shift; i <= text.size() - k; i ++) 
    {
        // std::cout << "i: " << i << "\n";

        unsigned kmer_new = encode_dna5(text.substr(i,k));

        if ((kmer_new & mask) > 0) {shift++; continue;} //kmer contains $

        // std::cout << "kmer_new: " << decode_dna5(kmer_new,k) << std::endl;
        unsigned interval_index = hashMap[kmer_new].lcp_interval_index;
        // printMap(k);
        lcp_interval temp_interval;

        if (lcpIntervals[interval_index] == std::nullopt) // Interval already computed or no valid kmer in Interval
        {
            // std::cout << "skipped 1" << std::endl;
            shift++;
            continue;
        }

        temp_interval = lcpIntervals[interval_index].value();
        unsigned count = rankSupport(temp_interval.right + 1) - rankSupport(temp_interval.left);

        
        unsigned long x = 0;
        unsigned long y = 0;
        
        // std::cout << "count: " << count << ", occurences: " << occurences << ", kmer: " <<text.substr(text_pos_kmer_new,k) <<  "\n";

        // std::cout << "text_pos_kmer_new: " << text_pos_kmer_new << ", text_pos_kmer_old + shift: " << text_pos_kmer_old + shift  << "\n";

        while (count > occurences && occurences > 0)
        {
            auto text_pos_kmer_new = suffixArray[temp_interval.left + x];
            auto text_pos_kmer_old = suffixArray[interval.left + y];
            

            if (text_pos_kmer_new != (text_pos_kmer_old + shift) )
            {   
                compressedSA.push_back(text_pos_kmer_new); 
                setValue(kmer_new, compressedSA.size() - 1, 1);
                x++;
                count--;
                std::cout << "x :" << x << ", y: " << y << std::endl;
                std::cout <<" interval.left + y: " << interval.left + y << " SuffixArray.size(): " << suffixArray.size() << std::endl;
                
                
                continue;
            }
            
            x++;
            y++;
            std::cout << "x :" << x << ", y: " << y << std::endl;
            std::cout << "Pattern for compression: " << text.substr(pat_pos_index,text.size()) << "\n";
            std::cout << "text_pos_kmerr_new: " << text.substr(text_pos_kmer_new,k) <<", " << text_pos_kmer_new << "\n";
            std::cout << "text_pos_kmerr_old + shift: " << text.substr(text_pos_kmer_old + shift,k) <<", " << text_pos_kmer_old + shift <<", shift: " << shift <<  "\n";
            std::cout <<" temp_interval.left + x: " << temp_interval.left + x << " SuffixArray.size(): " << suffixArray.size() << std::endl;
            std::cout <<" interval.left + y: " << interval.left + y << " SuffixArray.size(): " << suffixArray.size() << std::endl;
            std::cout << "FLAG" << std::endl;
            count--;
            occurences--;
        }

        if(occurences < count)
        {   
            while(count > 0)
            {   
                auto text_pos_kmer_new = suffixArray[temp_interval.left + x];
                compressedSA.push_back(text_pos_kmer_new);
                setValue(kmer_new, compressedSA.size() - 1, 1);
                x++;
                count--;
            }
            
        }

        // this -> printCompressedSA();
        setBits(temp_interval, 0); // Markiere alle Suffixe im Intervall als computed
        lcpIntervals[interval_index] = std::nullopt; // Setze das Intervall auf uninitialisiert
        setReferenceValue(kmer_new, shift,hashMap[kmer_old].occurences, hashMap[kmer_old].cSAindex, true);
        text_index ++;
        shift++;
    
    }
    

}

void CompressedSA::printIntervals(unsigned k)
{   
    unsigned long i = 0;
    for (const auto& interval : lcpIntervals)
    {
        if (!interval.has_value()) {
            std::cout<< "Index: " << i << ", Interval is uninitialized.\n";
            i++;
            continue; // Skip uninitialized intervals
        }
        std::cout <<"Index: " << i << ", Interval: " << "[" << interval->left
                  << ", " << interval->right << "],min: " << interval->min_index << ", priority: " << interval->priority << ", kmer: " << text.substr(suffixArray[interval->min_index], k) << "\n";
        i++;
    }
}

unsigned CompressedSA::getInterval(unsigned index)
{
    for (unsigned i = 0; i < lcpIntervals.size(); i++)
    {
        if (index >= lcpIntervals[i] -> left && index <= lcpIntervals[i] -> right)
        {
            return i;
        }
    }
    throw std::out_of_range("Kein Intervall gefunden");
}

unsigned CompressedSA::calc_priority(lcp_interval& interval) const
{   
    uint32_t quotient_lcp = 0;
    for (size_t i = interval.left + 1; i <= interval.right; i++) // left + 1 weil interval mit 4 Suffixen drei LCP_werte hat
    {
        quotient_lcp += lcpArray[i];
        std::cout << "lcpVal: " << lcpArray[i] << "\n";
    }
    quotient_lcp = std::ceil(quotient_lcp / (interval.right - interval.left + 1));
    return quotient_lcp;
}

unsigned CompressedSA::calc_priority2(lcp_interval& interval) const
{   
    std::cout << "FLAG" << std::endl;
    uint32_t quotient_lcp = 0;
    uint32_t minimum = std::numeric_limits<uint32_t>::max();
    std::vector<uint32_t> calculations;
    for (size_t i = interval.left + 1; i <= interval.right; i++) // left + 1 weil interval mit 4 Suffixen drei LCP_werte hat
    {
        if(lcpArray[i] < minimum) minimum = lcpArray[i];
        calculations.emplace_back((lcpArray[i]));
    }

    std::transform(calculations.begin(), calculations.end(), calculations.begin(), [minimum](uint32_t x) { return x - minimum; });
    quotient_lcp = std::ceil(std::accumulate(calculations.begin(), calculations.end(), 0) / calculations.size());
    std::cout << "quotient_lcp: " << quotient_lcp << "\n";
    return quotient_lcp;
 
}

void CompressedSA::initLCPintervalsAndHashmap(const unsigned k)
{
    unsigned i = sdsl::util::next_bit(computeSuffix,0);
    while (i < lcpArray.size())
    {
        lcp_interval interval = get_lcp_interval(i,k);
        interval.priority = calc_priority(interval);

        if (rankSupport(interval.right + 1) - rankSupport(interval.left) == 0)
        {
            // lcpIntervals.push_back(std::nullopt;);
            i = interval.right + 1;
            continue; // Skip if no valid k-mer in interval
        }
        unsigned long pos = suffixArray[interval.min_index];
        std::string kmer = text.substr(pos, k);
        uint64_t encodedKmer = encode_dna5(kmer);
        lcpIntervals.push_back(interval);
        this -> initHash(encodedKmer, lcpIntervals.size() - 1); // Store the interval index in the hash map
        
        i = interval.right + 1; // Move to the next interval

        for (unsigned long j = 0; j < kmer.size(); j++)
        {
            std::cout << kmer[j];
        }
        std::cout << "\n";
    }
}

void CompressedSA::runCompression(const unsigned k)
{   
    std::vector<unsigned> interval_indeces (lcpIntervals.size());
    std::iota(interval_indeces.begin(), interval_indeces.end(), 0);
    std::sort(interval_indeces.begin(), interval_indeces.end(), [this](unsigned i1, unsigned i2) { return lcpIntervals[i1] -> priority < lcpIntervals[i2] -> priority; });

    for (unsigned long i = 0; i < interval_indeces.size(); i++)
    {
        lcp_interval interval = lcpIntervals[interval_indeces[i]].value();
        std::string kmer = text.substr(suffixArray[interval.min_index], k);
        std::cout << "Index: " << interval_indeces[i] << ", Priority: " << lcpIntervals[interval_indeces[i]] -> priority << "kmer: "<< kmer << "\n";
    }

    unsigned prio_index;
    lcp_interval curr_interval;

    unsigned empty = 0;

    for (int i = interval_indeces.size() - 1; i >= 0; i--)
    {
        // std::cout << "interval.size(): " << interval_indeces.size() << "\n";
        if(lcpIntervals[interval_indeces[i]] == std::nullopt) 
        {
            interval_indeces.pop_back();    
            empty++; 
            continue;
        }

        prio_index = interval_indeces.back();
        interval_indeces.pop_back();
        curr_interval = lcpIntervals[prio_index].value();
        compression(k, curr_interval);

    }
    // std::cout << "isEmpty: " << interval_indeces.empty() << ", Intervals skipped: " << empty << "\n";
}

std::vector<int> CompressedSA::findPattern(std::string& kmer, unsigned k)
{
    if(kmer.size() != k) throw std::invalid_argument("Pattern size must be equal to kmer size");

    std::string retString = "[";

    std::vector<int> positions;
    bool isReference;
    bool isinCSA;

    int encoded_kmer = encode_dna5(kmer);
    hashValue curr_value = hashMap[encoded_kmer];

    isReference = (curr_value.refOccurrences != 0);
    isinCSA = (curr_value.occurences != 0);
    std::cout << "isReference: " << isReference << ", isinCSA: " << isinCSA << "\n";
    if(isinCSA)
    {   
        unsigned occ = curr_value.occurences;
        unsigned csa_index = curr_value.cSAindex; 
        std::cout << "csa_index: " << csa_index << ", occ: " << occ << "\n";
        for (unsigned long i = csa_index; i < csa_index + occ; i++)
        {
            positions.emplace_back(compressedSA[i]);
            retString += std::to_string(compressedSA[i]) + ",";
        }
    }

    if (isReference)
    {
        unsigned trace = curr_value.traceback_key;
        unsigned long refOcc = curr_value.refOccurrences;
        int shift = curr_value.shift;
        // std::cout << "trace: " << trace << ", refOcc: " << refOcc << ", shift: " << shift << "\n";

        for (unsigned long i = trace; i < trace + refOcc; i++)
        {
            positions.emplace_back(compressedSA[i] + shift);
            retString += std::to_string(compressedSA[i] + shift) + ",";
        }
    }
    
    retString[retString.size() - 1] = ']';

    std::cout << "Pattern kommt " << positions.size() << " mal vor in der Text, An Positionen: " << retString << "\n";
    
    return positions;
}
    
