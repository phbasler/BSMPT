cmake_minimum_required(VERSION 3.23)
project(PackageTest CXX)

find_package(BSMPT CONFIG REQUIRED)
find_package(OpenMP CONFIG QUIET)

add_executable(example src/example.cpp)
target_link_libraries(example BSMPT::Utility BSMPT::BounceSolution BSMPT::cmaes OpenMP::OpenMP_CXX
                      BSMPT::Models BSMPT::Minimizer)

target_link_libraries(example bsmpt::bsmpt)
