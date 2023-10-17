cmake_minimum_required(VERSION 3.12...3.15 FATAL_ERROR)
project(StaticMemTemplate LANGUAGES CXX)
add_compile_options(-std=c++14)

# open all compiler warnings
set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wpedantic")

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_executable(test_basic ./test/test_basic.cpp)
add_executable(test_blas ./test/test_blas.cpp)
add_executable(test_eigen ./test/test_eigen.cpp)
add_executable(test_kde ./test/test_kde.cpp)
add_executable(test_matrixdet ./test/test_matrixdet.cpp)

find_package (Eigen3 3.3 REQUIRED NO_MODULE)
add_executable(test_kdemv ./test/test_kdemv.cpp)
target_link_libraries (test_kdemv Eigen3::Eigen)

enable_testing()
add_test(NAME test_basic COMMAND ./test_basic)
add_test(NAME test_blas COMMAND ./test_blas)
add_test(NAME test_eigen COMMAND ./test_eigen)
add_test(NAME test_kde COMMAND ./test_kde)
add_test(NAME test_matrixdet COMMAND ./test_matrixdet)
add_test(NAME test_kdemv COMMAND ./test_kdemv)

