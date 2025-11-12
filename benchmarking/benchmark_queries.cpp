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
#include <future>    // std::future, std::async, std::launch
#include <sstream>   // parse_k_list



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

// Parse k-list like: "1,2,4-6,10"
static std::vector<unsigned> parse_k_list(const std::string& s) {
    std::vector<unsigned> out;
    std::string token;
    std::stringstream ss(s);
    while (std::getline(ss, token, ',')) {
        // trim spaces
        token.erase(0, token.find_first_not_of(" \t"));
        if (token.empty()) continue;
        token.erase(token.find_last_not_of(" \t") + 1);
        auto dash = token.find('-');
        if (dash == std::string::npos) {
            unsigned v = static_cast<unsigned>(std::stoul(token));
            out.push_back(v);
        } else {
            std::string a_str = token.substr(0, dash);
            std::string b_str = token.substr(dash + 1);
            // trim
            a_str.erase(0, a_str.find_first_not_of(" \t"));
            a_str.erase(a_str.find_last_not_of(" \t") + 1);
            b_str.erase(0, b_str.find_first_not_of(" \t"));
            b_str.erase(b_str.find_last_not_of(" \t") + 1);
            if (a_str.empty() || b_str.empty())
                throw std::invalid_argument("empty bound in range: '" + token + "'");
            unsigned a = static_cast<unsigned>(std::stoul(a_str));
            unsigned b = static_cast<unsigned>(std::stoul(b_str));
            if (a > b) std::swap(a, b);
            for (unsigned v = a; v <= b; ++v) out.push_back(v);
        }
    }
    if (out.empty()) throw std::invalid_argument("empty k_list");
    return out;
}

// -----------------------------
// Core worker: run one k (returns LookupResult)
// - Loads CSA for that k
// - Samples queries, measures SA/CSA in *nanoseconds*
// -----------------------------
static LookupResult run_one_k(const std::string& dataset_name,
                              const SuffixArray& SA,
                              unsigned k,
                              unsigned reps,
                              size_t num_queries,
                              bool shuffle_q) {
    using clock_t = std::chrono::steady_clock;

    // Load CSA specific to this k
    const std::string csa_file = "indices/" + dataset_name + "_k" + std::to_string(k) + "_csa.bin";
    auto CSA = std::make_unique<compressedSA>(compressedSA::load(csa_file));

    // Generate queries for this k
    std::vector<std::string> queries = findRandQueries(SA.getText(), k, num_queries);
    if (shuffle_q) {
        std::mt19937_64 rng(123456789);
        std::shuffle(queries.begin(), queries.end(), rng);
    }

    // Pre-allocate per-query accumulators (nanoseconds, averaged across reps)
    std::vector<double> sa_times_ns(num_queries, 0.0);
    std::vector<double> csa_times_ns(num_queries, 0.0);
    std::vector<size_t> sa_occurs(num_queries, 0);
    std::vector<size_t> csa_occurs(num_queries, 0);

    // Warmup
    for (const auto& q : queries) {
        (void)SA.search(q);
        std::string qm = q;
        (void)CSA->findPattern(qm, k);
    }

    const int R = static_cast<int>(reps);

    auto dur_ns = [](clock_t::time_point a, clock_t::time_point b) {
        return std::chrono::duration<double, std::micro>(b - a).count();
    };

    for (int r = 0; r < R; ++r) {
        std::vector<size_t> order(num_queries);
        std::iota(order.begin(), order.end(), size_t{0});
        std::mt19937_64 rng(static_cast<uint64_t>(123456789) + static_cast<uint64_t>(r));
        std::shuffle(order.begin(), order.end(), rng);

        const bool sa_first = (r % 2 == 0);

        for (size_t idx = 0; idx < num_queries; ++idx) {
            size_t i = order[idx];
            const auto& q = queries[i];

            auto measure_sa = [&]() -> double {
                auto t0 = clock_t::now();
                auto res = SA.search(q);
                auto t1 = clock_t::now();
                if (r == 0) sa_occurs[i] = res.size();
                return dur_ns(t0, t1);
            };

            auto measure_csa = [&]() -> double {
                std::string q_mut = q;
                auto t0 = clock_t::now();
                auto res = CSA->findPattern(q_mut, k);
                auto t1 = clock_t::now();
                if (r == 0) csa_occurs[i] = res.size();
                return dur_ns(t0, t1);
            };

            double dt_sa = 0.0, dt_csa = 0.0;
            if (sa_first) {
                dt_sa  = measure_sa();
                dt_csa = measure_csa();
            } else {
                dt_csa = measure_csa();
                dt_sa  = measure_sa();
            }

            sa_times_ns[i]  += dt_sa  / R;
            csa_times_ns[i] += dt_csa / R;

            if (r == 0 && sa_occurs[i] != csa_occurs[i]) {
                std::cerr << "WARNING: k=" << k << " query " << i << " (" << q
                          << ") mismatch: SA=" << sa_occurs[i]
                          << " CSA=" << csa_occurs[i] << "\n";
            }
        }
    }

    // Build LookupResult with means & medians (ns)
    double mean_sa   = mean(sa_times_ns);
    double mean_csa  = mean(csa_times_ns);
    double median_sa = percentile(sa_times_ns, 0.5);
    double median_csa= percentile(csa_times_ns, 0.5);

    return LookupResult(dataset_name, static_cast<uint32_t>(k),
                        mean_sa, mean_csa, median_sa, median_csa);
}

// -----------------------------
// main: multi-k parallel wrapper
// Usage:
//   prog <dataset_name> <k_list> <reps> [num_queries=100] [--shuffle]
// Example:
//   prog Human "8,10,12-16" 3 200 --shuffle
// -----------------------------
int main(int argc, char** argv) {
    bool shuffle_q = false;

    // Strip --shuffle flag
    {
        int w = 1;
        for (int r = 1; r < argc; ++r) {
            if (std::string(argv[r]) == "--shuffle") { shuffle_q = true; continue; }
            argv[w++] = argv[r];
        }
        argc = w;
    }

    if (argc < 4) {
        std::fprintf(stderr,
            "Usage: %s <dataset_name> <k_list> <reps> [num_queries=100] [--shuffle]\n",
            argv[0]);
        return 1;
    }

    const std::string dataset_name = argv[1];
    const std::string k_list_str   = argv[2];
    const unsigned reps            = static_cast<unsigned>(std::stoul(argv[3]));
    const size_t num_queries       = (argc >= 5) ? static_cast<size_t>(std::stoul(argv[4])) : 100;

    // Parse k list
    std::vector<unsigned> ks;
    try {
        ks = parse_k_list(k_list_str);
    } catch (const std::exception& e) {
        std::fprintf(stderr, "k_list parse error: %s\n", e.what());
        return 1;
    }

    printSeparator();
    std::cout << "QUERY BENCHMARK SETUP (multi-k)\n";
    printSeparator();
    std::cout << "Dataset:      " << dataset_name << "\n";
    std::cout << "k values:     ";
    for (size_t i = 0; i < ks.size(); ++i) std::cout << (i ? "," : "") << ks[i];
    std::cout << "\n";
    std::cout << "Reps/query:   " << reps << "\n";
    std::cout << "Num queries:  " << num_queries << "\n";
    std::cout << "Shuffle:      " << (shuffle_q ? "yes" : "no") << "\n\n";

    // Load SA once (shared by all tasks)
    std::unique_ptr<SuffixArray> SA;
    try {
        std::string sa_file = "indices/" + dataset_name + "_sa.bin";
        std::cout << "Loading SA from " << sa_file << "...\n";
        SA = std::make_unique<SuffixArray>(SuffixArray::load(sa_file));
        std::cout << "  Entries: " << SA->getSuffixArray().size() << "\n\n";
    } catch (...) {
        std::cerr << "ERROR loading SA\n";
        return 1;
    }

    // Launch per-k tasks in parallel
    std::vector<std::future<LookupResult>> futs;
    futs.reserve(ks.size());
    for (unsigned k : ks) {
        futs.emplace_back(std::async(std::launch::async, [&, k]() {
            return run_one_k(dataset_name, *SA, k, reps, num_queries, shuffle_q);
        }));
    }

    // Collect results
    std::vector<LookupResult> results;
    results.reserve(ks.size());
    for (auto& f : futs) results.push_back(f.get());

    // Print summary
    printSeparator();
    std::cout << "MULTI-k SUMMARY (ns/query; means & medians)\n";
    printSeparator();
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "dataset,k,mean_sa_ns,mean_csa_ns,median_sa_ns,median_csa_ns,slowdown\n";
    for (const auto& r : results) {
        std::cout << r.toLatexRow() << "\n";
    }
    printSeparator();
    std::cout << "DONE\n";
    printSeparator();
    return 0;
}
