cmake_minimum_required(VERSION 3.6)
project(CPPEx3)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra -Wall -Wvla -pthread -std=c++11")

set(SOURCE_FILES Complex.cpp GenericMatrixDriver.cpp BonusParallelChecker.cpp)
add_executable(CPPEx3 ${SOURCE_FILES})