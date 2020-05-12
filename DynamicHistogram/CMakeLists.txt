cmake_minimum_required(VERSION 3.10)
project(DynamicHistogram)

include_directories(src)

add_subdirectory(src)
add_subdirectory(test)

set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "Suppressing benchmark's tests" FORCE)
add_subdirectory(lib/benchmark)
add_subdirectory(lib/googletest)
