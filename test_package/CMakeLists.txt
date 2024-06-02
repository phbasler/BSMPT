cmake_minimum_required(VERSION 3.23)
project(PackageTest CXX)

find_package(BSMPT CONFIG REQUIRED)

add_executable(example src/example.cpp)
target_link_libraries(example bsmpt::Utility)