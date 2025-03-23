cmake_minimum_required(VERSION 3.10)
project(mathematical-algorithms)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED true)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_compile_options(-Wall -Wextra -Wpedantic)

add_executable(mathematical-algorithms src/main.cpp)

