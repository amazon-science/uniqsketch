cmake_minimum_required(VERSION 3.14)
project(UniqSketch VERSION 2.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Release by default
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# Dependencies. Threads: the query engine parses reads on std::threads
# (QuerySketchUtil.hpp), which every target that includes it needs to link.
find_package(ZLIB REQUIRED)
find_package(OpenMP)
find_package(Threads REQUIRED)

# Common include directories
include_directories(${CMAKE_SOURCE_DIR}/vendor ${CMAKE_SOURCE_DIR}/include)

# --- Main binaries ---

add_executable(uniqsketch src/uniqsketch.cpp)
target_link_libraries(uniqsketch PRIVATE ZLIB::ZLIB Threads::Threads)

add_executable(querysketch src/querysketch.cpp)
target_link_libraries(querysketch PRIVATE ZLIB::ZLIB Threads::Threads)

add_executable(comparesketch src/comparesketch.cpp)
target_link_libraries(comparesketch PRIVATE ZLIB::ZLIB)

if(OpenMP_CXX_FOUND)
    target_link_libraries(uniqsketch PRIVATE OpenMP::OpenMP_CXX)
    target_link_libraries(querysketch PRIVATE OpenMP::OpenMP_CXX)
    target_link_libraries(comparesketch PRIVATE OpenMP::OpenMP_CXX)
endif()

# --- Tests ---

enable_testing()

add_executable(uniqsketchtest test/uniqsketchtest.cpp)
target_link_libraries(uniqsketchtest PRIVATE ZLIB::ZLIB)

add_executable(querysketchtest test/querysketchtest.cpp)
target_link_libraries(querysketchtest PRIVATE ZLIB::ZLIB Threads::Threads)

add_executable(bloomfiltertest test/bloomfiltertest.cpp)
# The concurrency test for insert_make_change uses std::thread
target_link_libraries(bloomfiltertest PRIVATE Threads::Threads)

add_executable(comparesketchtest test/comparesketchtest.cpp)
target_link_libraries(comparesketchtest PRIVATE ZLIB::ZLIB)

# Builds a within-group sketch in-process (parallel job lists), so it needs OpenMP.
add_executable(treesketchtest test/treesketchtest.cpp)
target_link_libraries(treesketchtest PRIVATE ZLIB::ZLIB Threads::Threads)
if(OpenMP_CXX_FOUND)
    target_link_libraries(treesketchtest PRIVATE OpenMP::OpenMP_CXX)
endif()

# Compiles every *Util header into one translation unit; guards against
# reintroducing file-scope definitions that would stop the headers composing.
add_executable(headerstest test/headerstest.cpp)
target_link_libraries(headerstest PRIVATE ZLIB::ZLIB Threads::Threads)

# End-to-end: builds an index from the synthetic genomes in test/e2e/ and queries
# them in-process. Runs the pipeline functions, which carry OpenMP pragmas.
add_executable(e2etest test/e2etest.cpp)
target_link_libraries(e2etest PRIVATE ZLIB::ZLIB Threads::Threads)
if(OpenMP_CXX_FOUND)
    target_link_libraries(e2etest PRIVATE OpenMP::OpenMP_CXX)
endif()

set_target_properties(uniqsketchtest querysketchtest bloomfiltertest comparesketchtest treesketchtest headerstest e2etest
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/test)

# Tests run from the test/ source directory so they can find test data files
add_test(NAME uniqsketchtest COMMAND uniqsketchtest WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/test)
add_test(NAME querysketchtest COMMAND querysketchtest WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/test)
add_test(NAME bloomfiltertest COMMAND bloomfiltertest WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/test)
add_test(NAME comparesketchtest COMMAND comparesketchtest WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/test)
add_test(NAME treesketchtest COMMAND treesketchtest WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/test)
add_test(NAME headerstest COMMAND headerstest WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/test)
add_test(NAME e2etest COMMAND e2etest WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/test)
