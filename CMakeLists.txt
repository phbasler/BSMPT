cmake_minimum_required(VERSION 3.1)
project(BSMPT VERSION 1.0.0 LANGUAGES C CXX)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/tools/cmake")

set(CMAKE_CXX_STANDARD 11)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "RELEASE")
endif(NOT CMAKE_BUILD_TYPE)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")

find_package(Eigen3 REQUIRED)
if(EIGEN3_FOUND AND NOT TARGET Eigen3::Eigen)
  add_library(Eigen3::Eigen INTERFACE IMPORTED)
  set_target_properties(Eigen3::Eigen PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${EIGEN3_INCLUDE_DIR}")
endif()

if(NOT CMAES_ROOT)
  set(CMAES_ROOT "${CMAKE_SOURCE_DIR}/tools")
endif()
find_package(CMAES REQUIRED)

find_package(GSL REQUIRED)
find_package(OpenMP)

add_subdirectory(src)
