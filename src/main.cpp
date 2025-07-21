#include <iostream>
#include "suffix_array.hpp"
#include "compressV2.hpp"

int main()
{
    CompressedSA csa = CompressedSA("ATAACCGA$ATGACCGA$ATAACCGA$CTAACCGA$ATAACCGA");
    csa.initHashMap(2);
    csa.setReferenceValue(3, 1, 0, true);
    csa.setValue(3, 1, 3);
    hashValue hv = csa.getHashMapValue("AT");
    std::cout << "HashMap Value for AT: " << hv.kSAindex << ", " << hv.occurences << "\n";
    csa.printMap();
    return 0;
}

   