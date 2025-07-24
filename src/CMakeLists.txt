add_library(compressionTools
    fastaParser.cpp
    compressV2.cpp
    suffix_array.cpp
)

target_include_directories(compressionTools PUBLIC ${CMAKE_SOURCE_DIR}/include)