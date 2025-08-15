
#include "suffix_array.hpp"
#include "compressV2.hpp"
#include "fastaParser.hpp"
#include "test.hpp"
#include "compressedSA.hpp"
#include <memory>


using namespace sdsl;





int main(int argc, char* argv[])
{

    std::string filepath = argv[1];
    std::string kmer_to_find = argv[2];
    unsigned k = std::stoul(argv[3]);
    const std::string fastaData = parseFasta(filepath);
    std::cout << "datasize" <<fastaData.size() << "\n";  
    
    compressedSA e_csa (fastaData, k);
    std::vector<int> result = e_csa.findPattern(kmer_to_find, k);
    testCorrectness(fastaData, kmer_to_find, result);
}

   