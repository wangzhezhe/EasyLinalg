cmake_minimum_required(VERSION 3.12...3.15 FATAL_ERROR)
project(StaticMemTemplate LANGUAGES CXX C)
add_compile_options(-std=c++11)


add_executable(test test.cpp)