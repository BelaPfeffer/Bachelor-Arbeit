#include "compressV2.hpp"
#include "suffix_array.hpp"
#include <iostream>
#include "compressedSA.hpp"
#include <numeric>
#include <algorithm>
#include <cmath>
#include <iterator>  

#include <sdsl/suffix_array_algorithm.hpp>

uint64_t buildDollarMask(unsigned k) {
    uint64_t mask = 0;
    for (unsigned i = 0; i < k; ++i) {
        mask |= (uint64_t(4) << (3 * i));
    }
    return mask;
}

static inline uint8_t dna5_code(char c) {
    switch (c) { case 'A': return 0; case 'C': return 1; case 'G': return 2; case 'T': return 3; case '$': return 4; }
    throw std::invalid_argument("Invalid character in DNA sequence");
}

static inline uint64_t roll_kmer(uint64_t kmer, uint8_t next_code, unsigned k) {
    const uint64_t keep_mask = (k == 21) ? ~0ULL : ((1ULL << (3ULL * (k - 1))) - 1ULL);
    return ((kmer & keep_mask) << 3) | next_code;
}

void computeSA::setBits(lcp_interval& iv, bool bitVal)
{
    std::fill(computeSuffix.begin() + iv.left, computeSuffix.begin() + iv.right + 1, bitVal);
}


void computeSA::printComputeSA()
{
    std::cout << "Compressed Suffix Array:\n";
    for (unsigned long i = 0; i < CSA.size(); i++) {
        std::cout << i << "\t" << CSA[i] << "\t" << this -> text.substr(CSA[i]) << "\n";
    }
    std::cout << "\n";
}



void computeSA::printMap(uint64_t k)
{
    for (const auto& [key, value] : this -> hashMap) {

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

void computeSA::initHash(const uint64_t& kmer, unsigned long lcp_interval_index)
{
    this ->hashMap[kmer] = hashValue();
    this ->hashMap[kmer].lcp_interval_index = lcp_interval_index;

}


uint64_t computeSA::encode_dna5(const std::string& kmer)
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

std::string computeSA::decode_dna5(uint64_t encoded, unsigned k)
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

hashValue computeSA::getHashMapValue(const std::string& kmer)
{
    uint64_t encodedKmer = encode_dna5(kmer);
    return this->hashMap[encodedKmer];
}

void computeSA::setReferenceValue(const uint64_t encodedKmer, const int shift, unsigned long refOcurrences, const uint64_t traceback_key, const bool processed)
{
    this -> hashMap[encodedKmer].setReferenceValue(shift, refOcurrences, traceback_key, processed);
}


void computeSA::setValue(const uint64_t encodedKmer, const unsigned long cSAindex, const unsigned long occurences) // Index im komprimierten Suffix Array
{
   this -> hashMap[encodedKmer].setValue(cSAindex, occurences);
}

void computeSA::printSuffixArray() const 
{
        std::cout << "Suffix Array:\n";
        std::cout << "Index\tPosText\tLCP\tbitvektor\tSuffix\n";
        for (unsigned long i = 0; i < suffixArray.size(); i++) {
            std::cout << i << "\t"<<suffixArray[i] << "\t" << lcpArray[i] << "\t" << computeSuffix[i] << "\t" << text.substr(suffixArray[i]) << "\n";
        }
}

lcp_interval computeSA::get_lcp_interval(unsigned i, unsigned k)
{
    const sdsl::lcp_bitcompressed<>& lcp = lcpArray;
    unsigned long n = lcp.size() + 1; // Because LCP array has size n - 1
    unsigned long left = i;
    unsigned long right = i;
    unsigned long min_index = i;
    
    // Expand to the left
    while (left > 0 && lcp[left] >= k) {
        --left;
    }
    
    // Expand to the right
    while (right < n - 2 && lcp[right + 1] >= k) {
        ++right;
    }

    
    if (left != right)
    {
      min_index = left + 1;  
      for (unsigned long j = left + 1; j <= right; j++) {
          if (lcp[j] < lcp[min_index]) {
              min_index = j;
          }
      }   
    }
    
    // std::cout << "minLCPIndex: " << min_index << ", minLCPValue: " << lcp[min_index] << ", left: " << left << ", right: " << right << "\n";
    
    return lcp_interval(left, right, min_index);
}

//wird true wenn wert valid
bool computeSA::filterSA(const unsigned k, unsigned i)
{
    unsigned long suffix = suffixArray[i];
    if (suffix + k > text.size()) return false;
    for (unsigned j = 0; j < k; ++j)
        if (text[suffix + j] == '$') return false;
    return true;
}


void computeSA::initComputeSuffix(unsigned k)
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

void computeSA::compression(const unsigned k, lcp_interval& interval)
{
    // Buckets are disjoint; if an interval is present, it has content.
    if (interval.right < interval.left) return;

    // 1) Materialize the base interval into the CSA once.
    std::copy(suffixArray.begin() + interval.left,
              suffixArray.begin() + interval.right + 1,
              std::back_inserter(CSA));

    // Base pattern info
    const unsigned pat_pos_index  = suffixArray[interval.min_index];
    const unsigned pattern_length = lcpArray[interval.min_index];
    const unsigned long base_occ_total = interval.right - interval.left + 1;

    // Base k-mer (first window at pattern start), built once
    uint64_t kmer = 0;
    for (unsigned t = 0; t < k; ++t) kmer = (kmer << 3) | dna5_code(text[pat_pos_index + t]);
    const uint64_t kmer_old = kmer;

    // Record value for the base k-mer (now materialized in CSA).
    const unsigned long base_csa_index = CSA.size() - base_occ_total;
    setValue(kmer_old, base_csa_index, base_occ_total);

    // Invalidate the base interval entry in lcpIntervals (if present).
    const auto base_it = hashMap.find(kmer_old);
    if (base_it != hashMap.end() && base_it->second.lcp_interval_index < lcpIntervals.size())
        lcpIntervals[base_it->second.lcp_interval_index] = std::nullopt;

    // Mark every suffix in the base SA range as computed (optional; rank not used now).
    setBits(interval, 0);

    // Nothing more to do if LCP pattern is only k long.
    if (pattern_length <= k) return;

    // Reusable constants
    const uint64_t mask = buildDollarMask(k);

    // 2) Slide the k-mer window inside the longer pattern.
    for (unsigned long i = pat_pos_index + 1; i <= pat_pos_index + pattern_length - k; ++i)
    {
        const uint64_t shift = i - pat_pos_index;

        // O(1) rolling update for the new k-mer
        kmer = roll_kmer(kmer, dna5_code(text[i + k - 1]), k);
        const uint64_t kmer_new = kmer;

        // Skip kmers that include a sentinel ($).
        if ((kmer_new & mask) != 0) continue;

        // Must exist in hash map and lcpIntervals.
        auto it_new = hashMap.find(kmer_new);
        if (it_new == hashMap.end()) continue;

        const unsigned interval_index = it_new->second.lcp_interval_index;
        if (interval_index >= lcpIntervals.size()) continue;
        if (!lcpIntervals[interval_index].has_value()) continue;

        // Avoid self-interval reprocessing (same interval as base).
        if (base_it != hashMap.end() && interval_index == base_it->second.lcp_interval_index) continue;

        // Temporary interval for kmer_new
        lcp_interval temp_interval = lcpIntervals[interval_index].value();

        // Buckets are disjoint: number of items in this new interval:
        unsigned count = static_cast<unsigned>(temp_interval.right - temp_interval.left + 1);
        if (count == 0) { lcpIntervals[interval_index] = std::nullopt; continue; }

        unsigned long remaining_base = base_occ_total;
        unsigned long x = 0; // offset in temp_interval
        unsigned long y = 0; // offset in base interval

        // Two-pointer sweep while there are more new suffixes than base can cover
        while (count > remaining_base && remaining_base > 0) {
            const auto text_pos_kmer_new = suffixArray[temp_interval.left + x];
            const auto text_pos_kmer_old = suffixArray[interval.left + y];

            if (text_pos_kmer_new != (text_pos_kmer_old + shift)) {
                CSA.push_back(text_pos_kmer_new);
                setValue(kmer_new, CSA.size() - 1, 1);
                ++x; --count;
                continue;
            }
            ++x; ++y; --count; --remaining_base;
        }

        // Materialize the remainder if base cannot cover it.
        if (remaining_base < count) {
            while (count > 0) {
                const auto text_pos_kmer_new = suffixArray[temp_interval.left + x];
                CSA.push_back(text_pos_kmer_new);
                setValue(kmer_new, CSA.size() - 1, 1);
                ++x; --count;
            }
        }

        // Mark this new interval fully computed and invalidate its lcp entry.
        setBits(temp_interval, 0);
        lcpIntervals[interval_index] = std::nullopt;

        // Store reference metadata for kmer_new to the base k-mer (if different).
        if (kmer_new != kmer_old && base_it != hashMap.end()) {
            const auto& base_meta = base_it->second;
            setReferenceValue(kmer_new, static_cast<int>(shift),
                              base_meta.occurences, base_meta.cSAindex, true);
        }
    }
}


    




void computeSA::printIntervals(unsigned k)
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

unsigned computeSA::getInterval(unsigned index)
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

unsigned computeSA::calc_priority(lcp_interval& interval) const
{   
    uint32_t quotient_lcp = 0;
    for (size_t i = interval.left + 1; i <= interval.right; i++) // left + 1 weil interval mit 4 Suffixen drei LCP_werte hat
    {
        quotient_lcp += lcpArray[i];
        // std::cout << "lcpVal: " << lcpArray[i] << "\n";
    }
    quotient_lcp = std::ceil(quotient_lcp / (interval.right - interval.left + 1));
    return quotient_lcp;
}

unsigned computeSA::calc_priority2(lcp_interval& interval) const
{   
    // std::cout << "FLAG" << std::endl;
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
    // std::cout << "quotient_lcp: " << quotient_lcp << "\n";
    return quotient_lcp;
 
}

void computeSA::initLCPintervalsAndHashmap(const unsigned k)
{
    unsigned i = sdsl::util::next_bit(computeSuffix, 0);
    while (i < lcpArray.size())
    {
        lcp_interval interval = get_lcp_interval(i, k);
        interval.priority = calc_priority(interval);

        // One-time cheap check: any alive bit in [left,right]?
        unsigned alive = 0;
        for (unsigned p = interval.left; p <= interval.right; ++p) alive += computeSuffix[p];
        if (alive == 0) { i = interval.right + 1; continue; }

        unsigned long pos = suffixArray[interval.min_index];
        // encode directly from text (no substr needed here either, but substr is OK once)
        uint64_t encodedKmer = 0;
        for (unsigned t = 0; t < k; ++t) encodedKmer = (encodedKmer << 3) | dna5_code(text[pos + t]);

        lcpIntervals.push_back(interval);
        this->initHash(encodedKmer, lcpIntervals.size() - 1);

        i = interval.right + 1;
    }
}

void computeSA::runCompression(const unsigned k)
{
    std::vector<unsigned> interval_indeces(lcpIntervals.size());
    std::iota(interval_indeces.begin(), interval_indeces.end(), 0);

    // Guard against nullopt in comparator
    std::sort(interval_indeces.begin(), interval_indeces.end(),
        [this](unsigned i1, unsigned i2) {
            const auto &A = lcpIntervals[i1], &B = lcpIntervals[i2];
            if (!A) return false;
            if (!B) return true;
            return A->priority < B->priority; // ascending priority (your original)
        });

    // Process in reverse priority (largest first), safe reverse idiom
    for (size_t idx = interval_indeces.size(); idx-- > 0; ) {
        unsigned prio_index = interval_indeces[idx];
        auto &opt = lcpIntervals[prio_index];
        if (!opt) continue;

        compression(k, *opt);   // pass by reference (no copy)
        opt = std::nullopt;     // invalidate here (and in compression for the base)
    }
}



    
size_t computeSA::memoryUsageBytes() const {
    size_t totalMemory = 0;
    
    // 1. Original text string
    totalMemory += text.capacity(); // Allocated capacity, not just size()
    
    // 2. SDSL compressed suffix array
    totalMemory += size_in_bytes(suffixArray);
    
    // 3. SDSL LCP array
    totalMemory += size_in_bytes(lcpArray);
    
    // 4. Bit vector for computeSuffix
    totalMemory += size_in_bytes(computeSuffix);
    
    // 5. Rank support structure
    totalMemory += size_in_bytes(rankSupport);
    
    // 6. LCP intervals vector
    totalMemory += lcpIntervals.capacity() * sizeof(std::optional<lcp_interval>);
    // Add memory for the actual lcp_interval objects inside optionals
    for (const auto& opt_interval : lcpIntervals) {
        if (opt_interval.has_value()) {
            totalMemory += sizeof(lcp_interval);
        }
    }
    
    // 7. Hash map (more precise calculation)
    totalMemory += hashMap.size() * (sizeof(uint64_t) + sizeof(hashValue));
    // Hash map overhead: bucket array + node overhead
    totalMemory += hashMap.bucket_count() * sizeof(void*); // bucket pointers
    totalMemory += hashMap.size() * (sizeof(void*) + alignof(std::max_align_t)); // node overhead estimate
    
    // 8. Compressed SA vector
    totalMemory += CSA.capacity() * sizeof(int);
    
    // 9. Object overhead (approximate)
    totalMemory += sizeof(computeSA);
    
    return totalMemory;
}

compressedSA computeSA::exportSA () const
{
   compressedSA e_csa (this -> hashMap, this -> CSA, this -> text);
    return e_csa;
}


