cmake_minimum_required(VERSION 3.12)
project(LaminateOptimiser)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set output directories
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Find required packages
find_package(Eigen3 3.4 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIRS})

# Add GTest
include(FetchContent)
FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG release-1.12.1
)
FetchContent_MakeAvailable(googletest)

# First build all source libraries
add_subdirectory(src)

# Now add the tests
enable_testing()
add_subdirectory(examples)
add_subdirectory(tests)
add_subdirectory(benchmarks)
