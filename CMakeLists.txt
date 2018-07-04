cmake_minimum_required(VERSION 3.1)
project(BSMPT VERSION 1.0.0 LANGUAGES CXX)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/tools/cmake")

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")
set(TOOLSDIR ${CMAKE_SOURCE_DIR}/tools)

find_package(Eigen3 REQUIRED)
find_package(CMAES REQUIRED)
find_package(GSL REQUIRED)


add_subdirectory(src)
