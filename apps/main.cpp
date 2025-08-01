
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
    // csa_bitcompressed<> compSA;
    // lcp_bitcompressed<> compLCP;
    // construct_im(compSA, fastaData, 1);
    // construct_im(compLCP, fastaData, 1);
    // std::cout << "compSA size: " << compSA.size() << std::endl;
    // std::cout << "compLCP size: " << compLCP.size() << std::endl;
    // for (unsigned long i = 0; i < compSA.size(); i++) {
    //     std::cout << "i: " << i << ", char: " << compSA[i] << ", lcp: " << compLCP[i] << std::endl;
    // }

    
    CompressedSA csa (fastaData, k);
    csa.printSuffixArray();
    
    lcp_interval interval = csa.get_lcp_interval(14,k);

    std::cout << "Interval: " << interval.left << " " << interval.right << " " << interval.min_index << "\n";

    csa.compression(k, interval);
    csa.printMap(k);

    // csa.printSuffixArray();
    // csa.printMap(k);

}

   