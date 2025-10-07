
#include <benchmark/benchmark.h>
#include <sys/resource.h>
#include "compressedSA.hpp"
#include "fastaParser.hpp"
#include "suffix_array.hpp"
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
static std::vector<std::string> g_test_files;
static std::vector<std::string> g_parsed_data;
std::string kmer {};

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
    // 1. One-time setup (not measured) - now just get reference to pre-parsed data
    long file_index = state.range(0);
    const std::string& fastaData = g_parsed_data[file_index]; // Reference to pre-parsed data

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
    // 1. One-time setup (not measured) - now just get reference to pre-parsed data
    long file_index = state.range(0);
    unsigned k = state.range(1);
    const std::string& fastaData = g_parsed_data[file_index]; // Reference to pre-parsed data
    
    compressedSA csa (fastaData,k);
    state.counters["Exact Memory (Byte)"] = csa.memoryUsageBytes();
    state.counters["csasize"] = csa.csasize();
    // Measure memory before the measurement loop
    
    // 2. The actual measurement loop (measures time)
    for (auto _ : state) {
        // This is the operation whose time we want to measure.
        std::vector<int> result = csa.findPattern(kmer, k);
        // Prevent the compiler from optimizing away the object creation.
        benchmark::DoNotOptimize(result);
    }
    
    // 3. Set counters (after the loop!)
    // We report the absolute peak memory after execution.
    // Optional: Also report the difference if that's interesting to you.
    
}


// Function to pre-parse all FASTA files
void parseAllFiles() {
    g_parsed_data.reserve(g_test_files.size());
    
    for (const auto& filepath : g_test_files) {
        std::cout << "Parsing file: " << filepath << std::endl;
        std::string parsedData = parseFasta(filepath);
        g_parsed_data.push_back(std::move(parsedData));
        std::cout << "Parsed " << g_parsed_data.back().length() << " characters" << std::endl;
    }
}

// The corrected main() function
int main(int argc, char** argv) {
    // We need to separate the arguments from Google Benchmark.
    unsigned k = std::stoi(argv[argc - 1]);
    kmer = generate_random_sequence(k); // Default k=8 if not provided
    for (int i = 1; i < argc - 1; ++i) {
        if (argv[i][0] != '-') {
            g_test_files.push_back(argv[i]);
            // "Remove" the argument so Google Benchmark doesn't see it.
            // This is a common trick: overwrite it with the last one and shorten the list.
            argv[i] = argv[argc - 1];
            argc--;
            i--; // Check the new argument at this position again
        }
    }
    
    if (g_test_files.empty()) {
        fprintf(stderr, "Error: No input files specified.\n");
        return 1;
    }
    
    // Parse all files once before running benchmarks
    std::cout << "Pre-parsing all FASTA files..." << std::endl;
    parseAllFiles();
    std::cout << "Finished parsing all files." << std::endl;
    
    // Dynamic registration for each found file
    for (int i = 0; i < g_parsed_data.size(); ++i) {
        // Register the uncompressed suffix array benchmark
        benchmark::RegisterBenchmark("BM_uncompressedSA", &BM_uncompressedSA)->Arg(i);
        benchmark::RegisterBenchmark("BM_compressedSA", &BM_compressedSA)->Args({i, k});
   
        
        // You can easily add more benchmarks here that use the same parsed data
        // benchmark::RegisterBenchmark("BM_anotherBenchmark", &BM_anotherBenchmark)->Arg(i);
    }
    
    // Initialize and run Google Benchmark
    ::benchmark::Initialize(&argc, argv);
    ::benchmark::RunSpecifiedBenchmarks();
    ::benchmark::Shutdown();
    
    return 0;
}