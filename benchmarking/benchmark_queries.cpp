#include <benchmark/benchmark.h>
#include "compressedSA.hpp"
#include "suffix_array.hpp"
#include <random>
#include <vector>
#include <fstream>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <memory>

// Global state for loaded indices - USE POINTERS!
struct GlobalState {
    std::unique_ptr<SuffixArray> SA;        // Changed to pointer
    std::unique_ptr<compressedSA> CSA;      // Changed to pointer
    std::vector<std::string> queries;
    unsigned k;
    std::string dataset_name;
    bool loaded = false;
} g_state;

// Sample random k-mers from text
std::vector<std::string> sampleQueries(const std::string& text, unsigned k, 
                                       size_t num_queries, uint64_t seed) {
    if (text.size() < k) {
        throw std::runtime_error("Text too short for k-mer length");
    }
    
    std::mt19937_64 rng(seed);
    std::uniform_int_distribution<size_t> dist(0, text.size() - k);
    
    std::vector<std::string> queries;
    queries.reserve(num_queries);
    
    for (size_t i = 0; i < num_queries; ++i) {
        size_t pos = dist(rng);
        queries.push_back(text.substr(pos, k));
    }
    
    return queries;
}

// Benchmark uncompressed SA
static void BM_SA_Lookup(benchmark::State& state) {
    size_t query_idx = state.range(0);
    const std::string& query = g_state.queries[query_idx];
    
    for (auto _ : state) {
        auto result = g_state.SA->search(query);  // Use -> instead of .
        benchmark::DoNotOptimize(result);
        benchmark::ClobberMemory();
    }
}

// Benchmark compressed SA
static void BM_CSA_Lookup(benchmark::State& state) {
    size_t query_idx = state.range(0);
    std::string query = g_state.queries[query_idx];
    
    for (auto _ : state) {
        auto result = g_state.CSA->findPattern(query, g_state.k);  // Use -> instead of .
        benchmark::DoNotOptimize(result);
        benchmark::ClobberMemory();
    }
}

// Statistics helpers
double mean(const std::vector<double>& v) {
    if (v.empty()) return 0.0;
    return std::accumulate(v.begin(), v.end(), 0.0) / v.size();
}

double stddev(const std::vector<double>& v) {
    if (v.empty()) return 0.0;
    double m = mean(v);
    double sq_sum = 0.0;
    for (double x : v) sq_sum += (x - m) * (x - m);
    return std::sqrt(sq_sum / v.size());
}

double percentile(std::vector<double> v, double p) {
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    size_t idx = static_cast<size_t>(p * v.size());
    return v[std::min(idx, v.size() - 1)];
}

void printSeparator() {
    std::cout << "======================================\n";
}

int main(int argc, char** argv) {
    if (argc < 3) {
        std::fprintf(stderr, "Usage: %s <dataset_name> <k> [num_queries=100]\n", argv[0]);
        std::fprintf(stderr, "Example: %s kestrel 6 100\n", argv[0]);
        return 1;
    }
    
    g_state.dataset_name = argv[1];
    g_state.k = std::stoi(argv[2]);
    size_t num_queries = (argc >= 4) ? std::stoi(argv[3]) : 100;
    
    printSeparator();
    std::cout << "QUERY BENCHMARK SETUP\n";
    printSeparator();
    std::cout << "Dataset:       " << g_state.dataset_name << "\n";
    std::cout << "k-mer length:  " << g_state.k << "\n";
    std::cout << "Num queries:   " << num_queries << "\n";
    std::cout << "Random seed:   123456789\n";
    printSeparator();
    std::cout << "\n";
    
    // Construct filenames
    std::string sa_file = "indices/" + g_state.dataset_name + "_k" + 
                          std::to_string(g_state.k) + "_sa.bin";
    std::string csa_file = "indices/" + g_state.dataset_name + "_k" + 
                           std::to_string(g_state.k) + "_csa.bin";
    
    // Load indices using make_unique
    std::cout << "Loading uncompressed SA from " << sa_file << "...\n";
    try {
        g_state.SA = std::make_unique<SuffixArray>(SuffixArray::load(sa_file));
        std::cout << "  Entries: " << g_state.SA->getSuffixArray().size() << "\n";
        std::cout << "  Memory:  " << g_state.SA->memoryUsageBytes() / (1024.0 * 1024.0) << " MB\n";
    } catch (const std::exception& e) {
        std::cerr << "ERROR: Could not load SA: " << e.what() << "\n";
        std::cerr << "Did you run './construct_indices' first?\n";
        return 1;
    }
    
    std::cout << "\nLoading compressed SA from " << csa_file << "...\n";
    try {
        g_state.CSA = std::make_unique<compressedSA>(compressedSA::load(csa_file));
        std::cout << "  Entries: " << g_state.CSA->csasize() << "\n";
        std::cout << "  Memory:  " << g_state.CSA->memoryUsageBytes() / (1024.0 * 1024.0) << " MB\n";
    } catch (const std::exception& e) {
        std::cerr << "ERROR: Could not load CSA: " << e.what() << "\n";
        std::cerr << "Did you run './construct_indices' first?\n";
        return 1;
    }
    
    // Sample queries
    std::cout << "\nSampling " << num_queries << " random k-mers...\n";
    g_state.queries = sampleQueries(g_state.SA->getText(), g_state.k, num_queries, 123456789);
    std::cout << "  First query: " << g_state.queries[0] << "\n";
    std::cout << "  Last query:  " << g_state.queries[num_queries - 1] << "\n";
    
    g_state.loaded = true;
    
    std::cout << "\n";
    printSeparator();
    std::cout << "RUNNING GOOGLE BENCHMARK\n";
    printSeparator();
    std::cout << "\n";
    
    // Register benchmarks for each query
    for (size_t i = 0; i < num_queries; ++i) {
        benchmark::RegisterBenchmark("SA", BM_SA_Lookup)
            ->Args({static_cast<int64_t>(i)})
            ->UseRealTime()
            ->Unit(benchmark::kMicrosecond)
            ->Repetitions(5)
            ->ComputeStatistics("min", [](const std::vector<double>& v) -> double {
                return *std::min_element(v.begin(), v.end());
            })
            ->ComputeStatistics("max", [](const std::vector<double>& v) -> double {
                return *std::max_element(v.begin(), v.end());
            });
        
        benchmark::RegisterBenchmark("CSA", BM_CSA_Lookup)
            ->Args({static_cast<int64_t>(i)})
            ->UseRealTime()
            ->Unit(benchmark::kMicrosecond)
            ->Repetitions(5)
            ->ComputeStatistics("min", [](const std::vector<double>& v) -> double {
                return *std::min_element(v.begin(), v.end());
            })
            ->ComputeStatistics("max", [](const std::vector<double>& v) -> double {
                return *std::max_element(v.begin(), v.end());
            });
    }
    
    // Run Google Benchmark
    ::benchmark::Initialize(&argc, argv);
    ::benchmark::RunSpecifiedBenchmarks();
    ::benchmark::Shutdown();
    
    std::cout << "\n";
    printSeparator();
    std::cout << "POST-PROCESSING STATISTICS\n";
    printSeparator();
    std::cout << "\n";
    
    // Collect detailed measurements for statistics
    std::vector<double> sa_times, csa_times;
    std::vector<size_t> occurrences;
    
    std::cout << "Running single-pass measurement for statistics...\n";
    for (size_t i = 0; i < num_queries; ++i) {
        // Measure SA
        auto start = std::chrono::high_resolution_clock::now();
        auto sa_result = g_state.SA->search(g_state.queries[i]);
        auto end = std::chrono::high_resolution_clock::now();
        double sa_time = std::chrono::duration<double, std::micro>(end - start).count();
        sa_times.push_back(sa_time);
        
        // Measure CSA
        std::string query_copy = g_state.queries[i];
        start = std::chrono::high_resolution_clock::now();
        auto csa_result = g_state.CSA->findPattern(query_copy, g_state.k);
        end = std::chrono::high_resolution_clock::now();
        double csa_time = std::chrono::duration<double, std::micro>(end - start).count();
        csa_times.push_back(csa_time);
        
        occurrences.push_back(sa_result.size());
        
        // Verify correctness
        if (sa_result.size() != csa_result.size()) {
            std::cerr << "WARNING: Query " << i << " (" << g_state.queries[i] 
                      << ") has mismatched results!\n";
            std::cerr << "  SA occurrences:  " << sa_result.size() << "\n";
            std::cerr << "  CSA occurrences: " << csa_result.size() << "\n";
        }
    }
    
    std::cout << "Done.\n\n";
    
    // Print statistics
    std::cout << std::fixed << std::setprecision(3);
    
    printSeparator();
    std::cout << "UNCOMPRESSED SUFFIX ARRAY (SA)\n";
    printSeparator();
    std::cout << "Mean:      " << mean(sa_times) << " μs\n";
    std::cout << "Std Dev:   " << stddev(sa_times) << " μs\n";
    std::cout << "Median:    " << percentile(sa_times, 0.5) << " μs\n";
    std::cout << "P95:       " << percentile(sa_times, 0.95) << " μs\n";
    std::cout << "P99:       " << percentile(sa_times, 0.99) << " μs\n";
    std::cout << "Min:       " << *std::min_element(sa_times.begin(), sa_times.end()) << " μs\n";
    std::cout << "Max:       " << *std::max_element(sa_times.begin(), sa_times.end()) << " μs\n";
    
    std::cout << "\n";
    printSeparator();
    std::cout << "COMPRESSED SUFFIX ARRAY (CSA)\n";
    printSeparator();
    std::cout << "Mean:      " << mean(csa_times) << " μs\n";
    std::cout << "Std Dev:   " << stddev(csa_times) << " μs\n";
    std::cout << "Median:    " << percentile(csa_times, 0.5) << " μs\n";
    std::cout << "P95:       " << percentile(csa_times, 0.95) << " μs\n";
    std::cout << "P99:       " << percentile(csa_times, 0.99) << " μs\n";
    std::cout << "Min:       " << *std::min_element(csa_times.begin(), csa_times.end()) << " μs\n";
    std::cout << "Max:       " << *std::max_element(csa_times.begin(), csa_times.end()) << " μs\n";
    
    std::cout << "\n";
    printSeparator();
    std::cout << "PERFORMANCE COMPARISON\n";
    printSeparator();
    std::cout << "Slowdown factor (mean):   " << mean(csa_times) / mean(sa_times) << "x\n";
    std::cout << "Slowdown factor (median): " << percentile(csa_times, 0.5) / percentile(sa_times, 0.5) << "x\n";
    
    std::cout << "\n";
    
    // Save detailed CSV
    std::string csv_file = "results_" + g_state.dataset_name + "_k" + 
                           std::to_string(g_state.k) + ".csv";
    std::ofstream csv(csv_file);
    csv << "query_idx,query,sa_time_us,csa_time_us,occurrences\n";
    for (size_t i = 0; i < num_queries; ++i) {
        csv << i << "," << g_state.queries[i] << ","
            << sa_times[i] << "," << csa_times[i] << ","
            << occurrences[i] << "\n";
    }
    csv.close();
    
    std::cout << "Detailed results saved to: " << csv_file << "\n\n";
    
    // Output LaTeX table rows
    std::cout << "LaTeX table rows (for easy copy-paste):\n\n";
    std::cout << "% Lookup performance row:\n";
    std::cout << g_state.dataset_name << " & " << g_state.k << " & ... & ... & ... & ... & "
              << mean(sa_times) << " & " << mean(csa_times) << " \\\\\n\n";
    
    printSeparator();
    std::cout << "BENCHMARK COMPLETE!\n";
    printSeparator();
    
    return 0;
}