set(TEST ${CMAKE_PROJECT_NAME}_test)
file(GLOB_RECURSE TEST_SOURCES LIST_DIRECTORIES false *.h DynamicHistogramTest.cpp DynamicHistogram.cpp main.cpp)
add_executable(${TEST} ${TEST_SOURCES})
target_link_libraries(${TEST} PUBLIC ${CMAKE_PROJECT_NAME}_lib gmock gtest)
add_test(NAME ${TEST} COMMAND ${TEST})

set(BENCHMARK ${CMAKE_PROJECT_NAME}_benchmark)
file(GLOB_RECURSE BENCHMARK_SOURCES LIST_DIRECTORIES false *.h DynamicHistogramBenchmark.cpp)
add_executable(${BENCHMARK} ${BENCHMARK_SOURCES})
target_link_libraries(${BENCHMARK} PUBLIC ${CMAKE_PROJECT_NAME}_lib benchmark)

set(SOURCES ${TEST_SOURCES})

set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address -gmlt")
set (CMAKE_LINKER_FLAGS_DEBUG "${CMAKE_LINKER_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address -gmlt")

set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fno-omit-frame-pointer -gmlt")
set (CMAKE_LINKER_FLAGS_DEBUG "${CMAKE_LINKER_FLAGS_RELEASE} -fno-omit-frame-pointer -gmlt")
