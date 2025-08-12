#include "fastaParser.hpp"


std::string parseFasta (std::string path)
{
    gzFile file = gzopen(path.c_str(), "r");
    if (!file) {
        std::cerr << "Failed to open file\n";
        return {};
    }

    char buffer[1024];
    std::string current_id;
    std::string current_seq;
    std::vector<std::pair<std::string, std::string>> sequences;

    while (gzgets(file, buffer, sizeof(buffer))) {
        std::string line(buffer);

        if (line.back() == '\n') line.pop_back();

        if (line.empty()) continue;

        if (line[0] == '>') {
            if (!current_id.empty()) {
                sequences.emplace_back(current_id, current_seq);
                current_seq.clear();
            }
            current_id = line.substr(1); // skip '>'
        } else {
            current_seq += line;
        }
    }

    // Add the last sequence
    if (!current_id.empty()) {
        sequences.emplace_back(current_id, current_seq);
    }

    gzclose(file);

    // Print loaded sequences
    // for (auto &[id, seq] : sequences) {
    //     std::cout << ">" << id << "\n" << seq << "\n";
    // }

    std::string concatenatedText;
    for (const auto& [id, seq] : sequences) {
        concatenatedText += seq + "$"; // Append a special end character
    }
    // std::cout << "Concatenated text: " <<  concatenatedText << "\n";
    std::cout << "Parsing done. Concatenated text length: " << concatenatedText.length() << "\n";

    return concatenatedText;
}