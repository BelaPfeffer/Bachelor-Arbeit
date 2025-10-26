#include "compressedSA.hpp"
#include "suffix_array.hpp"
#include "fastaParser.hpp"
#include <iostream>
#include <filesystem>
#include <chrono>
#include <iomanip>

void printSeparator() {
    std::cout << "======================================\n";
}

int main(int argc, char** argv) {
    if (argc < 3) {
        std::fprintf(stderr, "Usage: %s <fasta_file> <k>\n", argv[0]);
        std::fprintf(stderr, "Example: %s data/kestrel.fa 6\n", argv[0]);
        return 1;
    }
    
    const std::string filepath = argv[1];
    const unsigned k = static_cast<unsigned>(std::stoi(argv[2]));
    
    // Extract dataset name from filepath
    std::filesystem::path p(filepath);
    std::string dataset_name = p.stem().string();
    
    // Create output directory
    std::filesystem::create_directories("indices");
    
    printSeparator();
    std::cout << "CONSTRUCTION PHASE\n";
    printSeparator();
    std::cout << "Dataset:  " << dataset_name << "\n";
    std::cout << "File:     " << filepath << "\n";
    std::cout << "k-mer:    " << k << "\n";
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
    
    // ===== Build Uncompressed Suffix Array =====
    printSeparator();
    std::cout << "BUILDING UNCOMPRESSED SUFFIX ARRAY\n";
    printSeparator();
    
    auto sa_start = std::chrono::high_resolution_clock::now();
    SuffixArray SA(text);
    auto sa_end = std::chrono::high_resolution_clock::now();
    auto sa_duration = std::chrono::duration_cast<std::chrono::milliseconds>(sa_end - sa_start);
    
    size_t sa_bytes = SA.memoryUsageBytes();
    size_t sa_entries = SA.getSuffixArray().size();
    
    std::cout << "Construction time:  " << sa_duration.count() / 1000.0 << " seconds\n";
    std::cout << "Memory usage:       " << sa_bytes << " bytes\n";
    std::cout << "Memory usage:       " << sa_bytes / (1024.0 * 1024.0) << " MB\n";
    std::cout << "SA entries:         " << sa_entries << "\n";
    std::cout << "Bytes per entry:    " << (double)sa_bytes / sa_entries << "\n";
    
    // Save to disk
    std::string sa_filename = "indices/" + dataset_name + "_k" + std::to_string(k) + "_sa.bin";
    std::cout << "\nSaving to: " << sa_filename << "\n";
    SA.save(sa_filename);
    
    size_t sa_filesize = std::filesystem::file_size(sa_filename);
    std::cout << "File size: " << sa_filesize / (1024.0 * 1024.0) << " MB\n\n";
    
    // ===== Build Compressed Suffix Array =====
    printSeparator();
    std::cout << "BUILDING COMPRESSED SUFFIX ARRAY\n";
    printSeparator();
    
    auto csa_start = std::chrono::high_resolution_clock::now();
    compressedSA CSA(text, k);
    auto csa_end = std::chrono::high_resolution_clock::now();
    auto csa_duration = std::chrono::duration_cast<std::chrono::milliseconds>(csa_end - csa_start);
    
    size_t csa_bytes = CSA.memoryUsageBytes();
    size_t csa_entries = CSA.csasize();
    
    std::cout << "Construction time:  " << csa_duration.count() / 1000.0 << " seconds\n";
    std::cout << "Memory usage:       " << csa_bytes << " bytes\n";
    std::cout << "Memory usage:       " << csa_bytes / (1024.0 * 1024.0) << " MB\n";
    std::cout << "CSA entries:        " << csa_entries << "\n";
    std::cout << "Bytes per entry:    " << (double)csa_bytes / csa_entries << "\n";
    
    // Save to disk
    std::string csa_filename = "indices/" + dataset_name + "_k" + std::to_string(k) + "_csa.bin";
    std::cout << "\nSaving to: " << csa_filename << "\n";
    CSA.save(csa_filename);
    
    size_t csa_filesize = std::filesystem::file_size(csa_filename);
    std::cout << "File size: " << csa_filesize / (1024.0 * 1024.0) << " MB\n\n";
    
    // ===== Summary =====
    printSeparator();
    std::cout << "SUMMARY\n";
    printSeparator();
    
    std::cout << std::fixed << std::setprecision(2);
    
    // Compression metrics
    double compression_ratio = (double)sa_entries / csa_entries;
    double memory_savings = (1.0 - (double)csa_bytes / sa_bytes) * 100.0;
    double construction_slowdown = (double)csa_duration.count() / sa_duration.count();
    
    std::cout << "Compression ratio:      " << compression_ratio << "x\n";
    std::cout << "Entry reduction:        " << sa_entries - csa_entries << " entries\n";
    std::cout << "Entry reduction:        " 
              << (1.0 - (double)csa_entries / sa_entries) * 100.0 << "%\n";
    std::cout << "Memory savings:         " << sa_bytes - csa_bytes << " bytes\n";
    std::cout << "Memory savings:         " << memory_savings << "%\n";
    std::cout << "Construction slowdown:  " << construction_slowdown << "x\n";
    
    printSeparator();
    std::cout << "\nCONSTRUCTION COMPLETE!\n\n";
    
    // Output LaTeX table row for easy copy-paste
    std::cout << "LaTeX table row:\n";
    std::cout << dataset_name << " & " << k << " & "
              << sa_bytes / (1024.0 * 1024.0) << " & "
              << csa_bytes / (1024.0 * 1024.0) << " & "
              << sa_entries / 1e6 << " & "
              << csa_entries / 1e6 << " \\\\\n";
    
    return 0;
}