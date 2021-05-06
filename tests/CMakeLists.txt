# SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas Müller
#
# SPDX-License-Identifier: GPL-3.0-or-later

add_library(Catch catch.cpp)
target_include_directories(Catch INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_features(Catch PUBLIC cxx_std_14)

set(TESTDIR ${CMAKE_CURRENT_SOURCE_DIR})

set(header
  GenerateTestCompares/C2HDM.h
  GenerateTestCompares/R2HDM.h
  GenerateTestCompares/RN2HDM.h
  )

set(src
  GenerateTestCompares/C2HDM.cpp
  GenerateTestCompares/R2HDM.cpp
  GenerateTestCompares/RN2HDM.cpp
  )

add_library(TestCompares ${header} ${src})
target_compile_features(TestCompares PUBLIC cxx_std_14)
target_link_libraries(TestCompares Minimizer)

file(GLOB SOURCE_FILES "Test-*.cpp" CONFIURE_DEPENDS)
add_executable(GenericTests ${SOURCE_FILES})
target_link_libraries(GenericTests Catch Minimizer Models Utility TestCompares)
target_compile_features(GenericTests PUBLIC cxx_std_14)

include(Catch)
catch_discover_tests(GenericTests)
