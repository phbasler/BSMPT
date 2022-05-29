# SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas Müller
#
# SPDX-License-Identifier: GPL-3.0-or-later

find_package(Catch2 3 QUIET)
if(NOT CATCH2_FOUND)
  include(FetchContent)
  FetchContent_Declare(
    Catch2
    GIT_SHALLOW    TRUE
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG        v3.0.1)
  FetchContent_MakeAvailable(Catch2)
endif()

list(APPEND CMAKE_MODULE_PATH ${Catch2_SOURCE_DIR}/extras)

add_subdirectory(catcheventlistener)
add_subdirectory(GenerateTestCompares)
add_subdirectory(unittests)
add_subdirectory(benchmarks)




