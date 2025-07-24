
#include "suffix_array.hpp"
#include "compressV2.hpp"
#include "fastaParser.hpp"


int main(int argc, char* argv[])
{
    std::string filepath = argv[1];
    unsigned k = std::stoul(argv[2]);
    std::string fastaData = parseFasta(filepath);
    CompressedSA csa (fastaData);
    csa.initHashMap(k);

    csa.printSuffixArray();
    csa.printMap(k);

}

   