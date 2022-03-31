# SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas Müller
#
# SPDX-License-Identifier: GPL-3.0-or-later

set(TESTDIR ${CMAKE_CURRENT_SOURCE_DIR})

set(SOURCE_FILES
  benchmark-ewpt-c2hdm.cpp)

if(TARGET Baryo)
  set(SOURCE_FILES
    ${SOURCE_FILES}
    benchmark-ewbg-c2hdm.cpp)
endif()
add_executable(benchmarks ${SOURCE_FILES})
target_link_libraries(benchmarks Catch Minimizer Models Utility TestCompares ThermalFunctions)
target_compile_features(benchmarks PUBLIC cxx_std_14)

if(TARGET Baryo)
  target_link_libraries(benchmarks Baryo)
endif()

include(Catch)
catch_discover_tests(benchmarks)
