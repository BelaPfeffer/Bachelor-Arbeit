#include <sdsl/suffix_arrays.hpp>
#include <sdsl/construct.hpp>
#include <sdsl/io.hpp>

#include "fastaParser.hpp"   // std::string parseFasta(path)

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

static inline void printSeparator() {
    std::cout << "======================================\n";
}

static std::string write_plain_text_for_sdsl(const std::string& dataset_name,
                                             const std::string& text) {
    std::filesystem::create_directories("indices");
    const std::string sdsl_in = "indices/" + dataset_name + "_sdsl_input.txt";
    std::ofstream out(sdsl_in, std::ios::binary);
    if (!out) throw std::runtime_error("cannot open output: " + sdsl_in);
    out.write(text.data(), static_cast<std::streamsize>(text.size()));
    out.close();
    return sdsl_in;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr, "Usage: %s <fasta_file>\n", argv[0]);
        return 1;
    }

    const std::string fasta_file = argv[1];
    const std::string dataset_name = std::filesystem::path(fasta_file).stem().string();
    const std::string out_idx = "indices/" + dataset_name + "_sdsl_csa_sada.sdsl";

    printSeparator();
    std::cout << "LOADING FASTA (for SDSL build)\n";
    printSeparator();

    auto tload0 = std::chrono::steady_clock::now();
    std::string text;
    try {
        text = parseFasta(fasta_file);
    } catch (const std::exception& e) {
        std::cerr << "ERROR parsing FASTA: " << e.what() << "\n";
        return 1;
    }
    auto tload1 = std::chrono::steady_clock::now();
    double load_s = std::chrono::duration<double>(tload1 - tload0).count();

    std::cout << "File:         " << fasta_file << "\n";
    std::cout << "Text length:  " << text.size() << " bp\n";
    std::cout << "Loading time: " << load_s << " s\n\n";

    printSeparator();
    std::cout << "SDSL CSA CONSTRUCTION\n";
    printSeparator();

    sdsl::csa_sada<> csa;
    auto t0 = std::chrono::steady_clock::now();
    try {
        const std::string sdsl_in = write_plain_text_for_sdsl(dataset_name, text);
        sdsl::construct(csa, sdsl_in, /*text_order=*/1);
    } catch (const std::exception& e) {
        std::cerr << "ERROR building SDSL CSA: " << e.what() << "\n";
        return 1;
    }
    auto t1 = std::chrono::steady_clock::now();
    double build_s = std::chrono::duration<double>(t1 - t0).count();

    sdsl::store_to_file(csa, out_idx);
    uint64_t mem_bytes  = sdsl::size_in_bytes(csa);
    uint64_t file_bytes = std::filesystem::file_size(out_idx);

    std::cout << "Construction time:  " << build_s << " seconds\n";
    std::cout << "Memory usage:       " << mem_bytes << " bytes ("
              << (mem_bytes / (1024.0*1024.0)) << " MB)\n";
    std::cout << "File saved:         " << out_idx << "\n";
    std::cout << "File size:          " << (file_bytes / (1024.0*1024.0)) << " MB\n";

    printSeparator();
    std::cout << "DONE (build)\n";
    printSeparator();
    return 0;
}
