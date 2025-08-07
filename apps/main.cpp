
#include "suffix_array.hpp"
#include "compressV2.hpp"
#include "fastaParser.hpp"


using namespace sdsl;


void checkCalcPos(const std::string& kmer, const std::vector<int>& positions, const std::string& text, const unsigned k)
{
    
    SuffixArray sa (text);
    std::vector<int> verified_pos = sa.search(kmer);
    bool found = false;

    for (int pos : positions)
    {
        found = (kmer == text.substr(pos, k));
    }

    if (found && verified_pos.size() == positions.size())
    {
        std::cout << "alle Pattern Richtig ermittlet " << "\n";
    }

    else if(found && verified_pos.size() != positions.size())
    {
        std::cout << "alle Pattern nicht richtig ermittelt aber zu wenig " << "\n";

    }
    else
    {
        std::cout << "Pattern nicht richtig ermittelt " << "\n";
    }
    return;
}


int main(int argc, char* argv[])
{

    std::string filepath = argv[1];
    std::string kmer_to_find = argv[2];
    unsigned k = std::stoul(argv[3]);
    const std::string fastaData = parseFasta(filepath);
    std::cout << "datasize" <<fastaData.size() << "\n";  

    

    
    CompressedSA csa (fastaData, k);
    
    
    csa.printSuffixArray();
    csa.printIntervals(k);
    csa.runCompression(k);
    csa.printMap(k);
    csa.printCompressedSA();

    std::vector<int> result = csa.findPattern(kmer_to_find, k);
    checkCalcPos(kmer_to_find, result, fastaData, k);

    
    // csa.printIntervals(k);
    
    // csa.printSuffixArray();
    

}

   