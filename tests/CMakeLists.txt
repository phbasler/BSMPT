# SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
# Müller
#
# SPDX-License-Identifier: GPL-3.0-or-later

find_package(Catch2 3 CONFIG QUIET)
find_package(GTest QUIET)

list(APPEND CMAKE_MODULE_PATH ${Catch2_SOURCE_DIR}/extras)

if(Catch2_FOUND)
  if(NOT "${Catch2_DIR}" STREQUAL "")
    set(CodeCoverageExcludesFromOtherPkgs ${CodeCoverageExcludesFromOtherPkgs}
                                          "${Catch2_DIR}/*")
  endif()
  if(NOT "${Catch2_ROOT_DIR}" STREQUAL "")
    set(CodeCoverageExcludesFromOtherPkgs ${CodeCoverageExcludesFromOtherPkgs}
                                          "${Catch2_ROOT_DIR}/*")
  endif()
endif()

add_subdirectory(catcheventlistener)
add_subdirectory(GenerateTestCompares)
add_subdirectory(unittests)
add_subdirectory(benchmarks)
