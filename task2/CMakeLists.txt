cmake_minimum_required(VERSION 3.21)
project(task2)

set(CMAKE_CXX_STANDARD 14)
find_package(OpenMP REQUIRED)
add_executable(task2 openmp.cpp)
target_link_libraries(task2 PRIVATE OpenMP::OpenMP_CXX)
