#include "compressedSA.hpp"
#include "suffix_array.hpp"
#include "test.hpp"      // for findRandQueries(...)
#include <random>
#include <vector>
#include <fstream>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <memory>
#include <chrono>
#include <iostream>

// -----------------------------
// Global state
// -----------------------------
struct GlobalState {
    std::unique_ptr<SuffixArray> SA;
    std::unique_ptr<compressedSA> CSA;
    std::vector<std::string> queries;
    unsigned k = 0;
    unsigned reps = 3;
    std::string dataset_name;
} g_state;

// -----------------------------
// Small helpers
// -----------------------------
static inline void printSeparator() {
    std::cout << "======================================\n";
}

static double mean(const std::vector<double>& v) {
    if (v.empty()) return 0.0;
    return std::accumulate(v.begin(), v.end(), 0.0) / v.size();
}

static double stddev(const std::vector<double>& v) {
    if (v.empty()) return 0.0;
    double m = mean(v);
    double sq = 0.0;
    for (double x : v) {
        double d = x - m;
        sq += d * d;
    }
    return std::sqrt(sq / v.size());
}

static double percentile(std::vector<double> v, double p) {
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    size_t idx = static_cast<size_t>(p * v.size());
    if (idx >= v.size()) idx = v.size() - 1;
    return v[idx];
}

int main(int argc, char** argv) {
    // Usage: prog <dataset_name> <k> <reps> [num_queries=100] [--shuffle]
    bool shuffle_q = false;

    // Strip shuffle flag if present
    {
        int w = 1;
        for (int r = 1; r < argc; ++r) {
            if (std::string(argv[r]) == "--shuffle") {
                shuffle_q = true;
                continue;
            }
            argv[w++] = argv[r];
        }
        argc = w;
    }

    if (argc < 4) {
        std::fprintf(stderr, "Usage: %s <dataset_name> <k> <reps> [num_queries=100] [--shuffle]\n", argv[0]);
        return 1;
    }

    g_state.dataset_name = argv[1];
    g_state.k = static_cast<unsigned>(std::stoul(argv[2]));
    g_state.reps = static_cast<unsigned>(std::stoul(argv[3]));
    size_t num_queries = (argc >= 5) ? static_cast<size_t>(std::stoul(argv[4])) : 100;

    printSeparator();
    std::cout << "QUERY BENCHMARK SETUP\n";
    printSeparator();
    std::cout << "Dataset:      " << g_state.dataset_name << "\n";
    std::cout << "k-mer length: " << g_state.k << "\n";
    std::cout << "Reps/query:   " << g_state.reps << "\n";
    std::cout << "Num queries:  " << num_queries << "\n";
    std::cout << "Random seed:  123456789\n";
    std::cout << "Shuffle:      " << (shuffle_q ? "yes" : "no") << "\n\n";

    // Filenames
    std::string sa_file  = "indices/" + g_state.dataset_name + "_k" + std::to_string(g_state.k) + "_sa.bin";
    std::string csa_file = "indices/" + g_state.dataset_name + "_k" + std::to_string(g_state.k) + "_csa.bin";

    // Load SA
    std::cout << "Loading SA from " << sa_file << "...\n";
    try {
        g_state.SA = std::make_unique<SuffixArray>(SuffixArray::load(sa_file));
        std::cout << "  Entries: " << g_state.SA->getSuffixArray().size() << "\n";
        std::cout << "  Memory:  " << g_state.SA->memoryUsageBytes() / (1024.0 * 1024.0) << " MB\n";
    } catch (...) {
        std::cerr << "ERROR loading SA\n";
        return 1;
    }

    // Load CSA
    std::cout << "Loading CSA from " << csa_file << "...\n";
    try {
        g_state.CSA = std::make_unique<compressedSA>(compressedSA::load(csa_file));
        std::cout << "  Entries: " << g_state.CSA->csasize() << "\n";
        std::cout << "  Memory:  " << g_state.CSA->memoryUsageBytes() / (1024.0 * 1024.0) << " MB\n";
    } catch (...) {
        std::cerr << "ERROR loading CSA\n";
        return 1;
    }

    // Generate queries
    std::cout << "\nSampling " << num_queries << " random k-mers...\n";
    g_state.queries = findRandQueries(g_state.SA->getText(), g_state.k, num_queries);

    if (shuffle_q) {
        std::mt19937_64 rng(123456789);
        std::shuffle(g_state.queries.begin(), g_state.queries.end(), rng);
    }

    using clock_t = std::chrono::steady_clock;

// Aggregatoren: direkt auf Größe setzen (nicht nur reserve),
// damit wir pro Query auf den Index akkumulieren können.
    std::vector<double> sa_times(num_queries, 0.0);
    std::vector<double> csa_times(num_queries, 0.0);
    std::vector<size_t> sa_occurs(num_queries, 0);
    std::vector<size_t> csa_occurs(num_queries, 0);

    printSeparator();
    std::cout << "MEASURING PER-QUERY TIMES (μs, averaged)\n";
    printSeparator();

    // Warmup (mutabler String für CSA, falls API non-const)
    for (const auto& q : g_state.queries) {
        (void)g_state.SA->search(q);
        std::string qm = q;
        (void)g_state.CSA->findPattern(qm, g_state.k);
    }

    const int R = static_cast<int>(g_state.reps);

    // Hilfsfunktion für gemessene Zeit (µs)
    auto dur_us = [](clock_t::time_point a, clock_t::time_point b){
        return std::chrono::duration<double, std::micro>(b - a).count();
    };

    // Drei (oder R) Pässe über die GESAMTE Query-Menge, jeweils anderer Shuffle.
    // Cross-over: in geraden Pässen SA→CSA, in ungeraden CSA→SA.
    for (int r = 0; r < R; ++r) {
        // deterministisches Shuffling pro Pass
        std::vector<size_t> order(num_queries);
        std::iota(order.begin(), order.end(), size_t{0});
        std::mt19937_64 rng(static_cast<uint64_t>(123456789) + static_cast<uint64_t>(r));
        std::shuffle(order.begin(), order.end(), rng);

        const bool sa_first = (r % 2 == 0);

        for (size_t idx = 0; idx < num_queries; ++idx) {
            size_t i = order[idx];
            const auto& q = g_state.queries[i];

            auto measure_sa = [&] () -> double {
                auto t0 = clock_t::now();
                auto res = g_state.SA->search(q);
                auto t1 = clock_t::now();
                if (r == 0) sa_occurs[i] = res.size(); // einmal pro Query erfassen
                return dur_us(t0, t1);
            };

            auto measure_csa = [&] () -> double {
                std::string q_mut = q; // falls API mutiert
                auto t0 = clock_t::now();
                auto res = g_state.CSA->findPattern(q_mut, g_state.k);
                auto t1 = clock_t::now();
                if (r == 0) csa_occurs[i] = res.size();
                return dur_us(t0, t1);
            };

            double dt_sa = 0.0, dt_csa = 0.0;
            if (sa_first) {
                dt_sa  = measure_sa();
                dt_csa = measure_csa();
            } else {
                dt_csa = measure_csa();
                dt_sa  = measure_sa();
            }

            // Über R-Pässe gemittelten Wert akkumulieren
            sa_times[i]  += dt_sa  / R;
            csa_times[i] += dt_csa / R;

            // Korrektheitscheck nur im ersten Pass (spart Arbeit)
            if (r == 0 && sa_occurs[i] != csa_occurs[i]) {
                std::cerr << "WARNING: Query " << i << " (" << q
                        << ") mismatch: SA=" << sa_occurs[i]
                        << " CSA=" << csa_occurs[i] << "\n";
            }
        }
    }


    // Summary
    printSeparator();
    std::cout << "SUMMARY (μs/query)\n";
    printSeparator();
    std::cout << std::fixed << std::setprecision(3);

    std::cout << "\nSA\n";
    std::cout << "Mean:   " << mean(sa_times) << " μs\n";
    std::cout << "StdDev: " << stddev(sa_times) << " μs\n";
    std::cout << "Median: " << percentile(sa_times, 0.5) << " μs\n";
    std::cout << "P95:    " << percentile(sa_times, 0.95) << " μs\n";
    std::cout << "P99:    " << percentile(sa_times, 0.99) << " μs\n";
    std::cout << "Min:    " << *std::min_element(sa_times.begin(), sa_times.end()) << " μs\n";
    std::cout << "Max:    " << *std::max_element(sa_times.begin(), sa_times.end()) << " μs\n";

    std::cout << "\nCSA\n";
    std::cout << "Mean:   " << mean(csa_times) << " μs\n";
    std::cout << "StdDev: " << stddev(csa_times) << " μs\n";
    std::cout << "Median: " << percentile(csa_times, 0.5) << " μs\n";
    std::cout << "P95:    " << percentile(csa_times, 0.95) << " μs\n";
    std::cout << "P99:    " << percentile(csa_times, 0.99) << " μs\n";
    std::cout << "Min:    " << *std::min_element(csa_times.begin(), csa_times.end()) << " μs\n";
    std::cout << "Max:    " << *std::max_element(csa_times.begin(), csa_times.end()) << " μs\n";

    std::cout << "\nSlowdown (mean): "
              << mean(csa_times) / mean(sa_times) << "x\n";

    // CSV export (write SA occurrences; CSA occurrences can be added if you want both)
    std::string csv_file = "results_" + g_state.dataset_name
                           + "_k" + std::to_string(g_state.k) + ".csv";
    std::ofstream csv(csv_file);
    csv << "query_idx,query,sa_time_us,csa_time_us,occurrences\n";
    for (size_t i = 0; i < num_queries; ++i) {
        csv << i << "," << g_state.queries[i] << ","
            << sa_times[i] << "," << csa_times[i] << ","
            << sa_occurs[i] << "\n";
    }
    csv.close();
    std::cout << "\nSaved: " << csv_file << "\n";

    printSeparator();
    std::cout << "DONE\n";
    printSeparator();

    return 0;
}
