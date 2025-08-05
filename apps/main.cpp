
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
    csa.printSuffixArray();
    


    // csa.findLCPintervals(k);


    
    lcp_interval interval = csa.get_lcp_interval(11,k);

    

    csa.compression(k, interval);
    csa.printMap(k);
    csa.printIntervals();
    // csa.printSuffixArray();
    

    // csa.printSuffixArray();
    // csa.printMap(k);

}

   