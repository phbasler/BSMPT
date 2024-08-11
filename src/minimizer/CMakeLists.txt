# SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
# Müller
#
# SPDX-License-Identifier: GPL-3.0-or-later

set(suffix "include/BSMPT/minimizer")
set(header_path "${BSMPT_SOURCE_DIR}/${suffix}")
set(header ${header_path}/Minimizer.h ${header_path}/MinimizePlane.h
           ${header_path}/MinimizeGSL.h)

add_library(Minimizer STATIC)
target_link_libraries(Minimizer PUBLIC Eigen3::Eigen GSL::gsl Threads::Threads
                                       Utility Models)
target_sources(
  Minimizer
  PUBLIC MinimizeGSL.cpp Minimizer.cpp MinimizePlane.cpp
  PUBLIC ${header})

if(cmaes_FOUND)
  add_subdirectory(LibCMAES)
  target_link_libraries(Minimizer PRIVATE Minimizer_CMAES)
endif()

if(NLopt_FOUND)
  add_subdirectory(LibNLOPT)
  target_link_libraries(Minimizer PUBLIC Minimizer_NLOPT)
endif()

target_include_directories(Minimizer PUBLIC ${BSMPT_SOURCE_DIR}/include
                                            ${BSMPT_BINARY_DIR}/include)
target_include_directories(Minimizer PUBLIC ${BSMPT_SOURCE_DIR}/include)
target_compile_features(Minimizer PUBLIC cxx_std_17)

set_property(TARGET Minimizer PROPERTY PUBLIC_HEADER ${header})

install(TARGETS Minimizer PUBLIC_HEADER DESTINATION ${suffix})
