cmake_minimum_required(VERSION 3.6)
project(GenMatDrive)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra -Wall -Wvla -pthread -std=c++11")

set(SOURCE_FILES Complex.cpp GenericMatrixDriver.cpp  Matrix.hpp)
add_executable(GenMatDrive ${SOURCE_FILES})

project(ParralelTest)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra -Wall -Wvla -pthread -O -DNDEBUG -std=c++11")

set(SOURCE_FILES Complex.cpp BonusParallelChecker.cpp  Matrix.hpp)
add_executable(ParTest ${SOURCE_FILES})
