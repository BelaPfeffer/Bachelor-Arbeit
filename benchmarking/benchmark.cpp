
#include <benchmark/benchmark.h>
#include <sys/resource.h>
#include "compressedSA.hpp"
#include "fastaParser.hpp"
#include "suffix_array.hpp"
#include "test.hpp"
#include <random>
#include <vector>
#include <string>

std::string generate_random_sequence(size_t k) {
    static const char nucleotides[] = {'A', 'C', 'G', 'T'};
    
    // Use a random device and Mersenne Twister engine
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dist(0, 3);

    std::string seq;
    seq.reserve(k);
    for (size_t i = 0; i < k; ++i) {
        seq.push_back(nucleotides[dist(gen)]);
    }
    return seq;
}



// Global variables for file paths and pre-parsed data
std::string text;
std::string kmer;
std::vector<uint64_t> result;

// Helper function for memory measurement
long getPeakMemoryUsageKB() {
    struct rusage usage;
    // RUSAGE_SELF queries memory for the current process
    if (getrusage(RUSAGE_SELF, &usage) == 0) {
        // On macOS ru_maxrss is in bytes, on Linux in kilobytes.
        // We need to consider this for portability.
#if defined(__APPLE__) && defined(__MACH__)
        return usage.ru_maxrss / 1024; // Convert from bytes to KB
#else
        return usage.ru_maxrss; // Already in KB
#endif
    }
    return 0;
}

// The corrected benchmark function
void BM_uncompressedSA(benchmark::State& state) {
    // std::cout << "Benchmarking uncompressedSA with kmer: " << kmer << std::endl;
    // 1. One-time setup (not measured) - now just get reference to pre-parsed data
    const std::string& fastaData = text;
    unsigned k = state.range(1); // Reference to pre-parsed data

    SuffixArray SA(fastaData);
    state.counters["Exact Memory (Byte)"] = SA.memoryUsageBytes();
    state.counters["sasize"] = SA.getSuffixArray().size();
    
    // Measure memory before the measurement loop

    
    // 2. The actual measurement loop (measures time)
    for (auto _ : state) {
        // This is the operation whose time we want to measure.
        std::vector<int> SA_result = SA.search(kmer);
        benchmark::DoNotOptimize(SA);
    }
   

    
   
    
}



void BM_compressedSA(benchmark::State& state) {
    // std::cout << "Benchmarking compressedSA with kmer: " << kmer << std::endl;
    // 1. One-time setup (not measured) - now just get reference to pre-parsed data
    unsigned k = state.range(0);
    const std::string& fastaData = text; // Reference to pre-parsed data
    
    compressedSA csa (fastaData,k);
    state.counters["Exact Memory (Byte)"] = csa.memoryUsageBytes();
    state.counters["csasize"] = csa.csasize();
    result = csa.findPattern(kmer, k);
    // Measure memory before the measurement loop
    
    // 2. The actual measurement loop (measures time)
    for (auto _ : state) {
        // This is the operation whose time we want to measure.
        std::vector<uint64_t> temp_result = csa.findPattern(kmer, k);
        // Prevent the compiler from optimizing away the object creation.
        benchmark::DoNotOptimize(temp_result);
    }
    
    // 3. Set counters (after the loop!)
    // We report the absolute peak memory after execution.
    // Optional: Also report the difference if that's interesting to you.
    
}


// Function to pre-parse all FASTA files
// void parseAllFiles() {
//     g_parsed_data.reserve(g_test_files.size());
    
//     for (const auto& filepath : g_test_files) {
//         std::cout << "Parsing file: " << filepath << std::endl;
//         std::string parsedData = parseFasta(filepath);
//         g_parsed_data.push_back(std::move(parsedData));
//         std::cout << "Parsed " << g_parsed_data.back().length() << " characters" << std::endl;
//     }
// }

// The corrected main() function
int main(int argc, char** argv) {
    // We need to separate the arguments from Google Benchmark.
    if (argc < 3) {
        fprintf(stderr, "Usage: %s <fasta_file1>  <k>\n", argv[0]);
        return 1;
    }
    unsigned k = std::stoi(argv[argc - 1]);
    std::string filepath = argv[1];
    text = parseFasta(filepath);
    kmer = findRandSequence(text,k);
    
    benchmark::RegisterBenchmark("BM_uncompressedSA", &BM_uncompressedSA)->Args({});
    benchmark::RegisterBenchmark("BM_compressedSA", &BM_compressedSA)->Args({k});
   
    
    // Initialize and run Google Benchmark
    ::benchmark::Initialize(&argc, argv);
    ::benchmark::RunSpecifiedBenchmarks();
    ::benchmark::Shutdown();
    std::cout << "===============================" << std::endl;
    std::cout << "BENCHMARKING DONE" << std::endl;
    std::cout << "===============================" << std::endl;
    std::cout << "Running correctness test..." << std::endl;
    testCorrectness(text, kmer, result);
    std::cout << "===============================" << std::endl;
    std::cout << "CORRECTNESS TEST DONE" << std::endl;
    std::cout << "===============================" << std::endl;
    
    return 0;
}