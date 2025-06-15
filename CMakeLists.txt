cmake_minimum_required(VERSION 3.15)
project(monte-carlo-options)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(-O3 -march=native -ffast-math)
endif()

add_executable(benchmarks pricing/benchmarks.cpp)
target_include_directories(benchmarks PRIVATE pricing)

find_package(Threads REQUIRED)
target_link_libraries(benchmarks Threads::Threads)
