
#include "suffix_array.hpp"
#include "compressV2.hpp"
#include "fastaParser.hpp"


using namespace sdsl;





int main(int argc, char* argv[])
{

    std::string filepath = argv[1];
    unsigned k = std::stoul(argv[2]);
    std::string fastaData = parseFasta(filepath);
    std::cout << "datasize" <<fastaData.size() << "\n";  

    
    CompressedSA csa (fastaData, k);
    csa.printMap(k);
    csa.printIntervals(k);
    // csa.printSuffixArray();
    csa.runCompression(k);
    // csa.printIntervals(k);
    
    // csa.printSuffixArray();
    

}

   