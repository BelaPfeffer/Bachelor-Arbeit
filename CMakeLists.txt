

cmake_minimum_required(VERSION 3.20...3.31)
project(seqan3_tutorial CXX)
set(CMAKE_CXX_STANDARD 17)

include(FetchContent)
FetchContent_Declare(
  seqan3
  GIT_REPOSITORY https://github.com/seqan/seqan3.git
  GIT_TAG        main  # or use a release tag like 3.3.0
)

FetchContent_MakeAvailable(seqan3)


FetchContent_Declare(
  sdsl
  GIT_REPOSITORY https://github.com/xxsds/sdsl-lite.git
  GIT_TAG        master
)

FetchContent_MakeAvailable(sdsl)

set(SDSL_INCLUDE_DIR ${sdsl_SOURCE_DIR}/include)
include_directories(${SDSL_INCLUDE_DIR})
# Pfad zu SeqAn3 (hier relativ zum Parent-Ordner, passe ggf. an)
# list(APPEND CMAKE_PREFIX_PATH "${CMAKE_CURRENT_SOURCE_DIR}/build/_deps/sdsl-src")

list(APPEND CMAKE_PREFIX_PATH "${CMAKE_CURRENT_SOURCE_DIR}/build/_deps/seqan3-src")
# SeqAn3 mit Versionsanforderung finden
find_package(seqan3 3.0 REQUIRED)



# Executable mit deinen Quellcodes
add_subdirectory(src)
add_subdirectory(benchmarking)
add_subdirectory(apps)
