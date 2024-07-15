cmake_minimum_required(VERSION 3.23)
project(PackageTest CXX)

find_package(BSMPT CONFIG REQUIRED)

add_compile_options(
  $<$<CXX_COMPILER_ID:MSVC>:/permissive->
  $<$<AND:$<CXX_COMPILER_ID:MSVC>,$<CONFIG:DEBUG>>:/bigobj>
  $<$<CXX_COMPILER_ID:MSVC>:/w44101>
  $<$<AND:$<CXX_COMPILER_ID:MSVC>,$<CONFIG:RELEASE>>:/Ox>)

add_executable(example src/example.cpp)
target_link_libraries(example BSMPT::Utility BSMPT::BounceSolution
                      BSMPT::Models BSMPT::Minimizer)
