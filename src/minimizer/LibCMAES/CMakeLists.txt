# SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
# Müller
#
# SPDX-License-Identifier: GPL-3.0-or-later

# include(CMAES)

set(suffix "include/BSMPT/minimizer/LibCMAES")
set(header_path "${BSMPT_SOURCE_DIR}/${suffix}")
set(header ${header_path}/MinimizeLibCMAES.h)

add_library(Minimizer_CMAES INTERFACE)
target_sources(Minimizer_CMAES INTERFACE MinimizeLibCMAES.cpp ${header})
target_link_libraries(Minimizer_CMAES INTERFACE Minimizer)

target_link_libraries(Minimizer_CMAES INTERFACE libcmaes::cmaes)

target_include_directories(
  Minimizer_CMAES INTERFACE ${BSMPT_SOURCE_DIR}/include
                            ${BSMPT_BINARY_DIR}/include)
target_include_directories(Minimizer_CMAES
                           INTERFACE ${BSMPT_SOURCE_DIR}/include)
target_compile_features(Minimizer_CMAES INTERFACE cxx_std_17)

set_property(TARGET Minimizer_CMAES PROPERTY PUBLIC_HEADER ${header})

install(TARGETS Minimizer_CMAES PUBLIC_HEADER DESTINATION ${suffix})
