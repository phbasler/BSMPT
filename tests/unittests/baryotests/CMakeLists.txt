# SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas Müller
#
# SPDX-License-Identifier: GPL-3.0-or-later

file(GLOB SOURCE_FILES "Test-*.cpp" CONFIURE_DEPENDS)
add_executable(BaryoTests ${SOURCE_FILES})
target_link_libraries(BaryoTests PRIVATE
  Catch2WithMain Minimizer Models Utility TestCompares ThermalFunctions Baryo)
target_compile_features(BaryoTests PUBLIC cxx_std_17)

include(Catch)
catch_discover_tests(BaryoTests)
