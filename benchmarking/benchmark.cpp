// #include <benchmark/benchmark.h>
// #include <sys/resource.h>
// #include "compressedSA.hpp"
// #include "fastaParser.hpp"
// #include "suffix_array.hpp"
// #include "test.hpp"
// #include <random>
// #include <vector>
// #include <string>
// #include <iostream>

// // Global inputs
// std::string text;
// std::string kmer;
// std::vector<uint64_t> result;

// // (Optional) peak RSS helper (currently unused)
// long getPeakMemoryUsageKB() {
//     struct rusage usage;
//     if (getrusage(RUSAGE_SELF, &usage) == 0) {
// #if defined(__APPLE__) && defined(__MACH__)
//         return usage.ru_maxrss / 1024; // bytes -> KB
// #else
//         return usage.ru_maxrss;        // already KB on Linux
// #endif
//     }
//     return 0;
// }

// static void BM_uncompressedSA(benchmark::State& state) {
//     const std::string& fastaData = text;

//     // Build once per benchmark invocation (not timed)
//     SuffixArray SA(fastaData);

//     // Report memory/size counters
//     state.counters["Exact Memory (Byte)"] = SA.memoryUsageBytes();
//     state.counters["sasize"] = SA.getSuffixArray().size();

//     for (auto _ : state) {
//         auto SA_result = SA.search(kmer);
//         benchmark::DoNotOptimize(SA_result);
//         benchmark::ClobberMemory();
//     }
// }

// static void BM_compressedSA(benchmark::State& state) {
//     unsigned k = static_cast<unsigned>(state.range(0));
//     const std::string& fastaData = text;

//     // Build once per benchmark invocation (not timed)
//     compressedSA csa(fastaData, k);

//     // Report memory/size counters
//     state.counters["Exact Memory (Byte)"] = csa.memoryUsageBytes();
//     state.counters["csasize"] = csa.csasize();

//     // Keep a reference result for correctness check after benchmarks
//     result = csa.findPattern(kmer, k);

//     for (auto _ : state) {
//         auto tmp = csa.findPattern(kmer, k);
//         benchmark::DoNotOptimize(tmp);
//         benchmark::ClobberMemory();
//     }
// }

// int main(int argc, char** argv) {
//     if (argc < 3) {
//         std::fprintf(stderr, "Usage: %s <fasta_file> <k>\n", argv[0]);
//         return 1;
//     }

//     const unsigned k = static_cast<unsigned>(std::stoi(argv[2]));
//     const std::string filepath = argv[1];

//     text = parseFasta(filepath);

//     // Reproducible query selection: choose one and print it
//     // (ensure findRandSequence uses a fixed or logged RNG seed)
//     kmer = findRandSequence(text, k);
//     std::cout << "Dataset: " << filepath << " (len=" << text.size()
//               << "), k=" << k << ", query=" << kmer << "\n";

//     benchmark::RegisterBenchmark("BM_uncompressedSA", &BM_uncompressedSA)
//         ->UseRealTime()
//         ->Unit(benchmark::kMicrosecond)
//         ->Repetitions(5)
//         ->ReportAggregatesOnly(true)
//         ->DisplayAggregatesOnly(true);

//     benchmark::RegisterBenchmark("BM_compressedSA", &BM_compressedSA)
//         ->Args({static_cast<int64_t>(k)})
//         ->UseRealTime()
//         ->Unit(benchmark::kMicrosecond)
//         ->Repetitions(5)
//         ->ReportAggregatesOnly(true)
//         ->DisplayAggregatesOnly(true);

//     ::benchmark::Initialize(&argc, argv);
//     ::benchmark::RunSpecifiedBenchmarks();
//     ::benchmark::Shutdown();

//     std::cout << "===============================\nBENCHMARKING DONE\n===============================\n";
//     std::cout << "Running correctness test...\n";
//     testCorrectness(text, kmer, result);
//     std::cout << "===============================\nCORRECTNESS TEST DONE\n===============================\n";
//     return 0;
// }
