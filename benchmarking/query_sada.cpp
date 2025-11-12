// bench_parallel.cpp
#include <sdsl/suffix_arrays.hpp>
#include <sdsl/io.hpp>
#include <sdsl/construct.hpp>

#include "fastaParser.hpp"   // std::string parseFasta(path)
#include "test.hpp"          // std::vector<std::string> findRandQueries(text, k, num)

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <cmath>

#include <thread>
#include <atomic>
#include <mutex>
#include <exception>

static inline void printSeparator() {
    std::cout << "======================================\n";
}
static double mean(const std::vector<double>& v) {
    if (v.empty()) return 0.0;
    return std::accumulate(v.begin(), v.end(), 0.0) / v.size();
}
static double variance(const std::vector<double>& v, double m) {
    if (v.size() < 2) return 0.0;
    double acc = 0.0;
    for (double x : v) {
        double d = x - m;
        acc += d * d;
    }
    return acc / (v.size() - 1);
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
        token.erase(0, token.find_first_not_of(" \t"));
        if (token.empty()) continue;
        token.erase(token.find_last_not_of(" \t") + 1);
        auto dash = token.find('-');
        if (dash == std::string::npos) {
            out.push_back(static_cast<unsigned>(std::stoul(token)));
        } else {
            std::string a = token.substr(0, dash), b = token.substr(dash + 1);
            a.erase(0, a.find_first_not_of(" \t"));
            a.erase(a.find_last_not_of(" \t") + 1);
            b.erase(0, b.find_first_not_of(" \t"));
            b.erase(b.find_last_not_of(" \t") + 1);
            unsigned av = static_cast<unsigned>(std::stoul(a));
            unsigned bv = static_cast<unsigned>(std::stoul(b));
            if (av > bv) std::swap(av, bv);
            for (unsigned v = av; v <= bv; ++v) out.push_back(v);
        }
    }
    if (out.empty()) throw std::invalid_argument("empty k_list");
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return out;
}

// LOCATE all occurrences (7-arg backward_search) and touch each SA entry
template<class CSA>
static inline size_t locate_all(const CSA& csa, const std::string& pat) {
    using size_type = typename CSA::size_type;
    if (pat.empty()) return 0;
    const size_type n = csa.size();
    if (n == 0) return 0;

    size_type l_in = 0, r_in = n - 1;
    size_type lb = 0, rb = 0;
    sdsl::backward_search(csa, l_in, r_in, pat.begin(), pat.end(), lb, rb);
    if (lb > rb || rb >= n) return 0;

    size_t occ = 0;
    for (size_type i = lb; i <= rb; ++i) {
        volatile auto pos = csa[i]; // emulate locate cost, prevent DCE
        (void)pos;
        ++occ;
    }
    return occ;
}

// Per-k benchmark (sequential; microseconds/query)
template<class CSA>
static std::tuple<double,double,double,double> bench_k(
    const CSA& csa,
    const std::string& text,
    unsigned k,
    unsigned reps,
    size_t num_queries,
    bool shuffle_q)
{
    using clock_t = std::chrono::steady_clock;

    auto queries = findRandQueries(text, k, num_queries);
    if (queries.size() < num_queries) {
        std::cerr << "WARNING: only " << queries.size() << " queries for k=" << k << "\n";
    }
    if (shuffle_q) {
        std::mt19937_64 rng(123456789);
        std::shuffle(queries.begin(), queries.end(), rng);
    }

    // warmup
    for (auto& q : queries) (void)locate_all(csa, q);

    std::vector<double> times_us(queries.size(), 0.0);
    const int R = static_cast<int>(reps);

    auto dur_us = [](clock_t::time_point a, clock_t::time_point b) {
        return std::chrono::duration<double, std::micro>(b - a).count();
    };

    for (int r = 0; r < R; ++r) {
        std::vector<size_t> order(queries.size());
        std::iota(order.begin(), order.end(), size_t{0});
        std::mt19937_64 rng(123456789ULL + r);
        std::shuffle(order.begin(), order.end(), rng);

        for (size_t j = 0; j < order.size(); ++j) {
            size_t i = order[j];
            auto t0 = clock_t::now();
            (void)locate_all(csa, queries[i]);
            auto t1 = clock_t::now();
            times_us[i] += dur_us(t0, t1) / R;
        }
    }

    double mean_us   = mean(times_us);
    double var_us    = variance(times_us, mean_us);
    double stddev_us = std::sqrt(var_us);
    double median_us = percentile(times_us, 0.5);
    return {mean_us, median_us, var_us, stddev_us};
}

int main(int argc, char** argv) {
    bool shuffle_q = false;
    int threads_arg = 0; // 0 = auto
    {
        int w = 1;
        for (int r = 1; r < argc; ++r) {
            std::string a = argv[r];
            if (a == "--shuffle") { shuffle_q = true; continue; }
            if (a == "--threads") {
                if (r + 1 >= argc) {
                    std::fprintf(stderr, "--threads requires a numeric argument\n");
                    return 1;
                }
                try {
                    threads_arg = std::stoi(argv[r + 1]);
                } catch (...) {
                    std::fprintf(stderr, "invalid --threads value: %s\n", argv[r + 1]);
                    return 1;
                }
                ++r; // skip the numeric argument
                continue;
            }
            argv[w++] = argv[r];
        }
        argc = w;
    }

    if (argc < 5) {
        std::fprintf(stderr,
            "Usage: %s <fasta_file> <k_list> <reps> <num_queries> [--shuffle] [--threads N]\n",
            argv[0]);
        return 1;
    }

    const std::string fasta_file = argv[1];
    const std::string k_list_str = argv[2];
    const unsigned reps          = static_cast<unsigned>(std::stoul(argv[3]));
    const size_t num_queries     = static_cast<size_t>(std::stoul(argv[4]));

    std::vector<unsigned> ks;
    try { ks = parse_k_list(k_list_str); }
    catch (const std::exception& e) {
        std::fprintf(stderr, "k_list parse error: %s\n", e.what());
        return 1;
    }

    const std::string dataset_name = std::filesystem::path(fasta_file).stem().string();
    const std::string index_path   = "indices/" + dataset_name + "_sdsl_csa_sada.sdsl";

    printSeparator();
    std::cout << "LOADING FASTA (for queries)\n";
    printSeparator();
    std::string text;
    try { text = parseFasta(fasta_file); }
    catch (const std::exception& e) {
        std::cerr << "ERROR parsing FASTA: " << e.what() << "\n";
        return 1;
    }
    std::cout << "File:        " << fasta_file << "\n";
    std::cout << "Text length: " << text.size() << " bp\n\n";

    printSeparator();
    std::cout << "LOADING SDSL INDEX\n";
    printSeparator();
    sdsl::csa_sada<> csa;
    if (!std::filesystem::exists(index_path)) {
        std::cerr << "ERROR: index not found at " << index_path << "\n";
        std::cerr << "Run the builder first: build_sdsl_csa <fasta_file>\n";
        return 1;
    }
    if (!sdsl::load_from_file(csa, index_path)) {
        std::cerr << "ERROR: failed to load index from " << index_path << "\n";
        return 1;
    }
    const uint64_t mem_bytes  = sdsl::size_in_bytes(csa);
    const uint64_t file_bytes = std::filesystem::file_size(index_path);
    std::cout << "Index loaded: " << index_path << "\n";
    std::cout << "Index bytes:  " << mem_bytes  << " (" << (mem_bytes  / (1024.0*1024.0)) << " MB)\n";
    std::cout << "File bytes:   " << file_bytes << " (" << (file_bytes / (1024.0*1024.0)) << " MB)\n\n";

    printSeparator();
    std::cout << "QUERY BENCHMARKS (LOCATE; microseconds per query)\n";
    printSeparator();
    std::cout << "Dataset:      " << dataset_name << "\n";
    std::cout << "k values:     ";
    for (size_t i = 0; i < ks.size(); ++i) std::cout << (i ? "," : "") << ks[i];
    std::cout << "\n";
    std::cout << "Reps/query:   " << reps << "\n";
    std::cout << "Num queries:  " << num_queries << "\n";
    std::cout << "Shuffle:      " << (shuffle_q ? "yes" : "no") << "\n";
    std::cout << "Threads arg:  " << (threads_arg > 0 ? std::to_string(threads_arg) : std::string("auto")) << "\n\n";

    struct Row { unsigned k; double mean_us; double median_us; double var_us; double stddev_us; };
    std::vector<Row> rows(ks.size());

    // --- parallel worker setup ------------------------------------------------
    std::atomic<size_t> next_idx{0};
    std::exception_ptr first_exc = nullptr;
    std::mutex exc_mu;

    auto worker = [&](void) {
        try {
            while (true) {
                size_t i = next_idx.fetch_add(1, std::memory_order_relaxed);
                if (i >= ks.size()) break;
                unsigned k = ks[i];

                // run benchmark for this k (reads shared csa)
                auto [mean_us, median_us, var_us, stddev_us] =
                    bench_k(csa, text, k, reps, num_queries, shuffle_q);

                rows[i] = {k, mean_us, median_us, var_us, stddev_us};
            }
        } catch (...) {
            std::lock_guard<std::mutex> lk(exc_mu);
            if (!first_exc) first_exc = std::current_exception();
        }
    };

    unsigned hw = std::thread::hardware_concurrency();
    size_t num_workers = 0;
    if (threads_arg > 0) {
        num_workers = static_cast<size_t>(threads_arg);
    } else if (hw > 0) {
        num_workers = static_cast<size_t>(hw);
    } else {
        num_workers = 4;
    }
    if (num_workers == 0) num_workers = 1;
    num_workers = std::min(num_workers, ks.size());

    std::vector<std::thread> threads;
    threads.reserve(num_workers);
    for (size_t t = 0; t < num_workers; ++t) threads.emplace_back(worker);

    for (auto& th : threads) th.join();

    if (first_exc) std::rethrow_exception(first_exc);

    // --- print results (single-threaded) -------------------------------------
    std::cout << "\n";
    printSeparator();
    std::cout << "SUMMARY\n";
    printSeparator();
    std::cout << "Index bytes (RAM): " << mem_bytes  << "  (" << (mem_bytes  / (1024.0*1024.0)) << " MB)\n";
    std::cout << "Index file bytes:  " << file_bytes << "  (" << (file_bytes / (1024.0*1024.0)) << " MB)\n\n";

    std::cout << "dataset,k,mean_csa_us,median_csa_us,var_csa_us,stddev_csa_us,index_bytes,index_file_bytes\n";
    std::cout << std::fixed << std::setprecision(2);
    for (auto& r : rows) {
        std::cout << dataset_name << "," << r.k << "," << r.mean_us << "," << r.median_us
                  << "," << r.var_us << "," << r.stddev_us
                  << "," << mem_bytes << "," << file_bytes << "\n";
    }

    printSeparator();
    std::cout << "DONE (bench)\n";
    printSeparator();
    return 0;
}
