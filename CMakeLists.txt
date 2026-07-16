cmake_minimum_required(VERSION 3.14)
project(UniqSketch VERSION 1.4.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Release by default
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# Dependencies
find_package(ZLIB REQUIRED)
find_package(OpenMP)

# Common include directories
include_directories(${CMAKE_SOURCE_DIR}/vendor ${CMAKE_SOURCE_DIR}/include)

# --- Main binaries ---

add_executable(uniqsketch src/uniqsketch.cpp)
target_link_libraries(uniqsketch PRIVATE ZLIB::ZLIB)

add_executable(querysketch src/querysketch.cpp)
target_link_libraries(querysketch PRIVATE ZLIB::ZLIB)

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
target_link_libraries(querysketchtest PRIVATE ZLIB::ZLIB)

add_executable(bloomfiltertest test/bloomfiltertest.cpp)

add_executable(comparesketchtest test/comparesketchtest.cpp)
target_link_libraries(comparesketchtest PRIVATE ZLIB::ZLIB)

set_target_properties(uniqsketchtest querysketchtest bloomfiltertest comparesketchtest
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/test)

# Tests run from the test/ source directory so they can find test data files
add_test(NAME uniqsketchtest COMMAND uniqsketchtest WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/test)
add_test(NAME querysketchtest COMMAND querysketchtest WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/test)
add_test(NAME bloomfiltertest COMMAND bloomfiltertest WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/test)
add_test(NAME comparesketchtest COMMAND comparesketchtest WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/test)
