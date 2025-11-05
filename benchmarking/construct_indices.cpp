#include "compressedSA.hpp"
#include "suffix_array.hpp"
#include "fastaParser.hpp"
#include "results.hpp"
#include <iostream>
#include <filesystem>
#include <chrono>
#include <iomanip>
#include <vector>
#include <thread>
#include <mutex>

std::mutex cout_mutex; // For thread-safe console output

void printSeparator() {
    std::cout << "======================================\n";
}

// Thread-safe logging
void log(const std::string& msg) {
    std::lock_guard<std::mutex> lock(cout_mutex);
    std::cout << msg;
}

// Worker function for each k-value
void buildCSAForK(const std::string& text, unsigned k, const std::string& dataset_name,
                  const MemoryResults& sa_result,
                  MemoryResults& csa_result,
                  std::chrono::milliseconds& duration_out) {
    
    log("[k=" + std::to_string(k) + "] Starting CSA construction...\n");
    
    auto csa_start = std::chrono::high_resolution_clock::now();
    compressedSA CSA(text, k);
    auto csa_end = std::chrono::high_resolution_clock::now();
    auto csa_duration = std::chrono::duration_cast<std::chrono::milliseconds>(csa_end - csa_start);
    
    // Set metadata
    csa_result.datasetName = dataset_name;
    csa_result.k = k;
    csa_result.numSAEntries = sa_result.numSAEntries;
    
    // Get CSA memory usage via setter method
    CSA.memoryUsageBytes(csa_result);
    
    // Save to disk
    std::string csa_filename = "indices/" + dataset_name + "_k" + std::to_string(k) + "_csa.bin";
    CSA.save(csa_filename);
    size_t csa_filesize = std::filesystem::file_size(csa_filename);
    
    // Calculate metrics
    double compression_ratio = (double)sa_result.numSAEntries / csa_result.numCSAEntries;
    double memory_savings = (1.0 - (double)csa_result.csaMemRaw / sa_result.saMemRaw) * 100.0;
    double entry_reduction = (1.0 - (double)csa_result.numCSAEntries / sa_result.numSAEntries) * 100.0;
    
    // Thread-safe logging
    {
        std::lock_guard<std::mutex> lock(cout_mutex);
        std::cout << "\n";
        printSeparator();
        std::cout << "CSA CONSTRUCTION COMPLETE (k=" << k << ")\n";
        printSeparator();
        std::cout << "Construction time:  " << csa_duration.count() / 1000.0 << " seconds\n";
        std::cout << "Memory usage:       " << csa_result.csaMemRaw << " bytes (" 
                  << csa_result.mbCSA() << " MB)\n";
        std::cout << "CSA entries:        " << csa_result.numCSAEntries << "\n";
        std::cout << "Bytes per entry:    " << (double)csa_result.csaMemRaw / csa_result.numCSAEntries << "\n";
        std::cout << "File saved:         " << csa_filename << "\n";
        std::cout << "File size:          " << csa_filesize / (1024.0 * 1024.0) << " MB\n";
        std::cout << "\nCompression metrics:\n";
        std::cout << "  Compression ratio:   " << std::fixed << std::setprecision(2) 
                  << compression_ratio << "x\n";
        std::cout << "  Entry reduction:     " << entry_reduction << "%\n";
        std::cout << "  Memory savings:      " << memory_savings << "%\n";
        std::cout << std::endl;
    }
    
    duration_out = csa_duration;
}

int main(int argc, char** argv) {
    if (argc < 3) {
        std::fprintf(stderr, "Usage: %s <fasta_file> <k1> [k2 k3 ...] [--threads N]\n", argv[0]);
        std::fprintf(stderr, "Example: %s data/kestrel.fa 6 8 10\n", argv[0]);
        std::fprintf(stderr, "         %s data/kestrel.fa 6 8 10 --threads 4\n", argv[0]);
        std::fprintf(stderr, "         %s data/kestrel.fa 6\n", argv[0]);
        return 1;
    }

    const std::string filepath = argv[1];
    
    // Parse all k values and thread count
    std::vector<unsigned> k_values;
    unsigned num_threads = std::thread::hardware_concurrency(); // Default to hardware concurrency
    
    for (int i = 2; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--threads" && i + 1 < argc) {
            num_threads = std::stoi(argv[++i]);
        } else {
            k_values.push_back(static_cast<unsigned>(std::stoi(argv[i])));
        }
    }
    
    // Extract dataset name from filepath
    std::filesystem::path p(filepath);
    std::string dataset_name = p.stem().string();
    
    // Create output directory
    std::filesystem::create_directories("indices");

    printSeparator();
    std::cout << "CONSTRUCTION PHASE (PARALLEL)\n";
    printSeparator();
    std::cout << "Dataset:  " << dataset_name << "\n";
    std::cout << "File:     " << filepath << "\n";
    std::cout << "k-values: ";
    for (size_t i = 0; i < k_values.size(); ++i) {
        std::cout << k_values[i];
        if (i < k_values.size() - 1) std::cout << ", ";
    }
    std::cout << "\n";
    std::cout << "Threads:  " << num_threads << " (max: " 
              << std::thread::hardware_concurrency() << ")\n";
    printSeparator();
    std::cout << "\n";
    
    // Load FASTA
    std::cout << "Loading FASTA file...\n";
    auto load_start = std::chrono::high_resolution_clock::now();
    std::string text = parseFasta(filepath);
    auto load_end = std::chrono::high_resolution_clock::now();
    auto load_time = std::chrono::duration_cast<std::chrono::seconds>(load_end - load_start);
    std::cout << "  Text length:   " << text.size() << " bp\n";
    std::cout << "  Loading time:  " << load_time.count() << " seconds\n\n";

    // ===== Build/Load Uncompressed Suffix Array (ONCE) =====
    printSeparator();
    std::cout << "SUFFIX ARRAY (SHARED FOR ALL k)\n";
    printSeparator();
    
    // Check if SA already exists for this dataset
    std::string sa_filename = "indices/" + dataset_name + "_sa.bin";
    SuffixArray SA(text);
    MemoryResults sa_result;
    sa_result.datasetName = dataset_name;
    auto sa_duration = std::chrono::milliseconds(0);
    
    if (std::filesystem::exists(sa_filename)) {
        std::cout << "Found existing SA file: " << sa_filename << "\n";
        std::cout << "Loading from disk...\n";
        auto load_sa_start = std::chrono::high_resolution_clock::now();
        SA.load(sa_filename);
        auto load_sa_end = std::chrono::high_resolution_clock::now();
        sa_duration = std::chrono::duration_cast<std::chrono::milliseconds>(load_sa_end - load_sa_start);
        std::cout << "Loading time:       " << sa_duration.count() / 1000.0 << " seconds\n";
    } else {
        std::cout << "Building new suffix array...\n";
        auto sa_start = std::chrono::high_resolution_clock::now();
        // SA is already constructed above with text
        auto sa_end = std::chrono::high_resolution_clock::now();
        sa_duration = std::chrono::duration_cast<std::chrono::milliseconds>(sa_end - sa_start);
        std::cout << "Construction time:  " << sa_duration.count() / 1000.0 << " seconds\n";
        
        // Save for future use
        std::cout << "Saving to: " << sa_filename << "\n";
        SA.save(sa_filename);
    }
    
    // Get SA memory usage via setter
    SA.memoryUsageBytes(sa_result);
    size_t sa_filesize = std::filesystem::file_size(sa_filename);
    
    std::cout << "Memory usage:       " << sa_result.saMemRaw << " bytes (" 
              << sa_result.saMemRaw / (1024.0 * 1024.0) << " MB)\n";
    std::cout << "SA entries:         " << sa_result.numSAEntries << "\n";
    std::cout << "Bytes per entry:    " << (double)sa_result.saMemRaw / sa_result.numSAEntries << "\n";
    std::cout << "File size:          " << sa_filesize / (1024.0 * 1024.0) << " MB\n\n";

    // ===== Build Compressed Suffix Arrays in Parallel =====
    printSeparator();
    std::cout << "BUILDING CSA STRUCTURES IN PARALLEL\n";
    printSeparator();
    std::cout << "Launching " << k_values.size() << " parallel construction(s)...\n\n";
    
    auto parallel_start = std::chrono::high_resolution_clock::now();
    
    // Storage for results - initialize each with a copy of sa_result
    std::vector<MemoryResults> csa_results(k_values.size(), sa_result);
    std::vector<std::chrono::milliseconds> csa_durations(k_values.size());
    
    // Launch threads in batches to respect thread limit
    std::vector<std::thread> threads;
    size_t k_idx = 0;
    
    while (k_idx < k_values.size()) {
        // Launch batch of threads
        size_t batch_size = std::min(num_threads, static_cast<unsigned>(k_values.size() - k_idx));
        threads.clear();
        
        for (size_t i = 0; i < batch_size; ++i) {
            size_t current_idx = k_idx + i;
            threads.emplace_back(buildCSAForK, 
                                std::cref(text), 
                                k_values[current_idx], 
                                std::cref(dataset_name),
                                std::cref(sa_result),
                                std::ref(csa_results[current_idx]),
                                std::ref(csa_durations[current_idx]));
        }
        
        // Wait for batch to complete
        for (auto& thread : threads) {
            thread.join();
        }
        
        k_idx += batch_size;
    }
    
    auto parallel_end = std::chrono::high_resolution_clock::now();
    auto total_parallel_time = std::chrono::duration_cast<std::chrono::milliseconds>(parallel_end - parallel_start);

    // ===== Summary Table =====
    std::cout << "\n";
    printSeparator();
    std::cout << "SUMMARY\n";
    printSeparator();
    std::cout << std::fixed << std::setprecision(2);
    
    std::cout << "\nDataset: " << dataset_name << "\n";
    std::cout << "Text size: " << text.size() << " bp\n";
    std::cout << "SA entries: " << sa_result.numSAEntries << "\n";
    std::cout << "SA memory: " << sa_result.saMemRaw / (1024.0 * 1024.0) << " MB\n";
    std::cout << "Total parallel construction time: " << total_parallel_time.count() / 1000.0 << " seconds\n\n";
    
    std::cout << std::left << std::setw(6) << "k" 
              << std::right << std::setw(12) << "CSA entries"
              << std::setw(12) << "CSA MB"
              << std::setw(12) << "Ratio"
              << std::setw(12) << "Reduction%"
              << std::setw(12) << "Time(s)\n";
    std::cout << std::string(66, '-') << "\n";
    
    for (size_t i = 0; i < k_values.size(); ++i) {
        unsigned k = k_values[i];
        const auto& csa_result = csa_results[i];
        auto csa_duration = csa_durations[i];
        
        double compression_ratio = (double)sa_result.numSAEntries / csa_result.numCSAEntries;
        double reduction_pct = (1.0 - (double)csa_result.numCSAEntries / sa_result.numSAEntries) * 100.0;
        
        std::cout << std::left << std::setw(6) << k
                  << std::right << std::setw(12) << csa_result.numCSAEntries
                  << std::setw(12) << csa_result.mbCSA()
                  << std::setw(12) << compression_ratio
                  << std::setw(12) << reduction_pct
                  << std::setw(12) << (csa_duration.count() / 1000.0) << "\n";
    }
    
    printSeparator();
    std::cout << "\nCONSTRUCTION COMPLETE!\n\n";
    
    // Output LaTeX table rows
    std::cout << "LaTeX table rows:\n";
    for (size_t i = 0; i < k_values.size(); ++i) {
        csa_results[i].printLatexRow();
    }

    return 0;
}