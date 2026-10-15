# SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
# Müller
#
# SPDX-License-Identifier: GPL-3.0-or-later

set(TESTDIR ${CMAKE_CURRENT_SOURCE_DIR})

file(GLOB SOURCE_FILES CONFIGURE_DEPENDS "Test-*.cpp")
add_executable(GenericTests ${SOURCE_FILES})
target_link_libraries(
  GenericTests
  PRIVATE GSL::gsl
          Catch2::Catch2WithMain
          CatchEventListener
          Minimizer
          MinimumTracer
          Models
          Utility
          TestCompares
          ThermalFunctions
          BounceSolution
          TransitionTracer
          GW
          BaryoFHCK
          VacuumProfile)

if(BSMPTBuildExecutables)
  target_sources(GenericTests PRIVATE CheckExecutables.cpp)
  target_compile_definitions(GenericTests
                             PRIVATE BASE_PATH="${CMAKE_SOURCE_DIR}")

endif()

target_compile_definitions(GenericTests
                           PRIVATE TEST_DATA_PATH="${CMAKE_CURRENT_SOURCE_DIR}")

target_compile_features(GenericTests PUBLIC cxx_std_17)

include(Catch)
catch_discover_tests(GenericTests)
