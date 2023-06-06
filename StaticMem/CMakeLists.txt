cmake_minimum_required(VERSION 3.12...3.15 FATAL_ERROR)
project(StaticMemTemplate LANGUAGES CXX C)
add_compile_options(-std=c++11)


add_executable(test_static_3by3 test_static_3by3.cpp)
add_executable(test_static_4by4 test_static_4by4.cpp)
add_executable(test_static_8by8 test_static_8by8.cpp)