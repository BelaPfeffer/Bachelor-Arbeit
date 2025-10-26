# --- Google Benchmark Integration ---
include(FetchContent)
FetchContent_Declare(
  benchmark
  GIT_REPOSITORY https://github.com/google/benchmark.git
  GIT_TAG        v1.8.3
)
set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "" FORCE)
set(BENCHMARK_ENABLE_INSTALL OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(benchmark)

# --- Original Benchmark Executable ---
add_executable(benchmark_script benchmark.cpp)
target_compile_definitions(benchmark_script PRIVATE BENCHMARK_MODE)
target_link_libraries(benchmark_script
  PRIVATE
    compressionTools
    seqan3::seqan3
    benchmark::benchmark
    benchmark::benchmark_main
)

# --- NEW: Construction Executable (no benchmark dependency) ---
add_executable(construct_indices construct_indices.cpp)
target_link_libraries(construct_indices
  PRIVATE
    compressionTools
    seqan3::seqan3
)

# --- NEW: Query Benchmark Executable ---
add_executable(benchmark_queries benchmark_queries.cpp)
target_link_libraries(benchmark_queries
  PRIVATE
    compressionTools
    seqan3::seqan3
    benchmark::benchmark
)