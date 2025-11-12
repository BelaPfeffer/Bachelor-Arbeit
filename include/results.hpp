#ifndef MEMORY_RESULTS_HPP
#define MEMORY_RESULTS_HPP
#pragma once
#include <sstream>
#include <string>
#include <cstdint>
#include <iostream>
#include <iomanip>


struct MemoryResults {
    // Metadata
    std::string datasetName;
    unsigned k;

    // Entry counts
    uint64_t numSAEntries;     // SA entries (usually = text length)
    uint64_t numCSAEntries;    // CSA entries (usually = number of k-mer intervals)

    // SA (Suffix Array) memory (in bytes)
    uint64_t saMemRaw;         // Raw SA memory (e.g., using uint64_t entries)
    uint64_t saMem4B;          // Theoretical SA memory at 4 bytes per entry (32 bit)
    uint64_t lcpMemRaw;        // Raw LCP array memory
    uint64_t saLcpMem;         // (SA 4B + LCP) – only if LCP is needed for queries

    // CSA + Hashmap memory (in bytes)
    uint64_t csaMemRaw;
    uint64_t csaMem4B;        // Theoretical CSA memory at 4 bytes per entry (32 bit)
    uint64_t hashmapMemRaw;    // Raw memory for hashmap
    uint64_t csaHMTotalMem;    // Total: CSA + Hashmap

    // Repetitiveness metrics
    double avgMinLCP = 0.0;              // Average min LCP per k-mer interval
    double pctIntervalsMinLCPgeK = 0.0;  // % of intervals with minLCP >= k

      // ---------- Default constructor ----------
    MemoryResults()
        : datasetName(""),
          k(0),
          numSAEntries(0),
          numCSAEntries(0),
          saMemRaw(0),
          saMem4B(0),
          lcpMemRaw(0),
          saLcpMem(0),
          csaMemRaw(0),
          hashmapMemRaw(0),
          csaHMTotalMem(0),
          avgMinLCP(0.0),
          pctIntervalsMinLCPgeK(0.0)
    {}

    // ---------- Copy constructor ----------
    MemoryResults(const MemoryResults& other)
        : datasetName(other.datasetName),
          k(other.k),
          numSAEntries(other.numSAEntries),
          numCSAEntries(other.numCSAEntries),
          saMemRaw(other.saMemRaw),
          saMem4B(other.saMem4B),
          lcpMemRaw(other.lcpMemRaw),
          saLcpMem(other.saLcpMem),
          csaMemRaw(other.csaMemRaw),
          hashmapMemRaw(other.hashmapMemRaw),
          csaHMTotalMem(other.csaHMTotalMem),
          avgMinLCP(other.avgMinLCP),
          pctIntervalsMinLCPgeK(other.pctIntervalsMinLCPgeK)
    {}

    
    // Setter for SA + LCP data
    void setSAData(uint64_t saRawBytes, uint64_t lcpRawBytes, uint64_t saEntries) {
        saMemRaw = saRawBytes;
        lcpMemRaw = lcpRawBytes;
        numSAEntries = saEntries;
        saMem4B = numSAEntries * 4;
        saLcpMem = saMem4B + lcpMemRaw;
    }

    // Setter for CSA + Hashmap data
    void setCSAData(uint64_t csaRawBytes, uint64_t hashmapRawBytes, uint64_t csaEntries) {
        csaMemRaw = csaRawBytes;
        csaMem4B = csaEntries * 4;
        hashmapMemRaw = hashmapRawBytes;
        numCSAEntries = csaEntries;
        csaHMTotalMem = csaMem4B + hashmapMemRaw;
    }

    // Setter for repetitiveness data
    void setRepetitivenessData(double avgMin, double pctMinGeK) {
        avgMinLCP = avgMin;
        pctIntervalsMinLCPgeK = pctMinGeK;
    }

    // Convenience methods for memory in MB
    double mbSA() const { return static_cast<double>(saMem4B) / 1e6; }
    double mbSAplusLCP() const { return static_cast<double>(saLcpMem) / 1e6; }
    double mbCSA() const { return static_cast<double>(csaMem4B) / 1e6; }
    double mbCSAplusHM() const { return static_cast<double>(csaHMTotalMem) / 1e6; }

    // Compression factors
    double compFactorCSA() const { return static_cast<double>(saMem4B) / csaMem4B; }
    double compFactorTotal() const { return static_cast<double>(saLcpMem) / csaHMTotalMem; }

    // Print a LaTeX-formatted table row
    void printLatexRow() const {
        std::cout << std::fixed << std::setprecision(2);
        std::cout << " & "
            << k << " & "
            << mbSA() << " & "
            << mbCSA() << " & "
            << (numSAEntries / 1e6) << " & "
            << (numCSAEntries / 1e6) << " & "
            << compFactorCSA() << " & "
            << mbSAplusLCP() << " & "
            << mbCSAplusHM() << " & "
            << compFactorTotal()
            << " \\\\" << std::endl;
    }
};



struct LookupResult {
    // === Metadata ===
    std::string dataset;  // e.g. "Human", "Cod"
    uint32_t k;           // k-mer size

    // === Lookup performance (in nanoseconds) ===
    double meanSaNs;      // Mean lookup time for SA
    double meanCsaNs;     // Mean lookup time for CSA + Hashmap
    double medianSaNs;    // (optional) Median lookup time for SA
    double medianCsaNs;   // (optional) Median lookup time for CSA

    // === Derived metrics ===
    double slowdown;      // CSA / SA (mean)
    double speedup;       // SA / CSA (if you want the inverse)

    // === Constructor ===
    LookupResult(std::string dataset_, uint32_t k_,
                 double meanSaNs_, double meanCsaNs_,
                 double medianSaNs_ = 0, double medianCsaNs_ = 0)
        : dataset(std::move(dataset_)), k(k_),
          meanSaNs(meanSaNs_), meanCsaNs(meanCsaNs_),
          medianSaNs(medianSaNs_), medianCsaNs(medianCsaNs_)
    {
        slowdown = (meanSaNs > 0) ? (meanCsaNs / meanSaNs) : 0.0;
        speedup = (meanCsaNs > 0) ? (meanSaNs / meanCsaNs) : 0.0;
    }

    // === Helper: print as LaTeX table row (optional) ===
    std::string toLatexRow() const {
        std::ostringstream ss;
        ss << std::fixed << std::setprecision(2);
        ss << " & " << k << " & "
           << meanSaNs << " & " << meanCsaNs << " & "
           << slowdown << " \\\\";
        return ss.str();
    }
};


#endif // MEMORY_RESULTS_HPP
