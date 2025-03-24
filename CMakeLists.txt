cmake_minimum_required(VERSION 3.10)
project(mathematical-algorithms)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED true)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -funroll-loops")

add_compile_options(-Wall -Wextra -Wpedantic -O2 -funroll-loops)

add_executable(karatsuba src/karatsuba.cpp)
add_executable(pseudorandom src/pseudorandom.cpp)
add_executable(types src/types.cpp)
