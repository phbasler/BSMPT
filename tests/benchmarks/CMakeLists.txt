# SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
# Müller
#
# SPDX-License-Identifier: GPL-3.0-or-later

find_package(benchmark QUIET)

if(benchmark_FOUND)
  set(TESTDIR ${CMAKE_CURRENT_SOURCE_DIR})

  set(SOURCE_FILES benchmark-ewpt-c2hdm.cpp)

  if(BSMPTCompileBaryo)
    set(SOURCE_FILES ${SOURCE_FILES} benchmark-ewbg-c2hdm.cpp)
  endif()
  add_executable(benchmarks ${SOURCE_FILES})
  target_link_libraries(benchmarks Minimizer Models Utility TestCompares
                        ThermalFunctions)
  target_link_libraries(benchmarks benchmark::benchmark)
  target_compile_features(benchmarks PUBLIC cxx_std_14)

  if(BSMPTCompileBaryo)
    target_link_libraries(benchmarks Baryo)
  endif(BSMPTCompileBaryo)

else()
  message(
    WARNING "Did not find google benchmark. Benchmarks will not be compiled.")
endif()
