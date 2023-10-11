cmake_minimum_required(VERSION 3.12...3.15 FATAL_ERROR)
project(StaticMemTemplate LANGUAGES CXX C)
add_compile_options(-std=c++11)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_executable(test_basic ./test/test_basic.cpp)
add_executable(test_blas ./test/test_blas.cpp)
add_executable(test_eigen ./test/test_eigen.cpp)
add_executable(test_kde ./test/test_kde.cpp)
add_executable(test_matrixdet ./test/test_matrixdet.cpp)

find_package (Eigen3 3.3 REQUIRED NO_MODULE)
add_executable(test_kdemv ./test/test_kdemv.cpp)
target_link_libraries (test_kdemv Eigen3::Eigen)
