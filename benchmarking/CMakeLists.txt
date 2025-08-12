
  # --- Google Benchmark Integration (THE FIX) ---
  # We use FetchContent to compile benchmark from source, ensuring it uses the
  # same compiler and standard library as the rest of our project.

  # This should ideally be in your top-level CMakeLists.txt, but it's fine here too.
  include(FetchContent)

  FetchContent_Declare(
    benchmark
    GIT_REPOSITORY https://github.com/google/benchmark.git
    GIT_TAG        v1.8.3 # Using a specific version is good practice
  )

  # Prevent gbenchmark from building its own tests or being installed globally
  set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "" FORCE)
  set(BENCHMARK_ENABLE_INSTALL OFF CACHE BOOL "" FORCE)

  # This makes the `benchmark` targets available to our project
  FetchContent_MakeAvailable(benchmark)

  # --- Your Benchmark Executable ---

  # This line remains the same
  add_executable(benchmark_script benchmark.cpp)

  target_compile_definitions(benchmark_script PRIVATE BENCHMARK_MODE)

  # This line also remains the same. The target names provided by FetchContent
  # are identical to the ones from find_package, making the switch seamless.
  target_link_libraries(benchmark_script
    PRIVATE
      compressionTools
      seqan3::seqan3
      benchmark::benchmark
      benchmark::benchmark_main
  )

