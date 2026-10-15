# SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
# Müller
#
# SPDX-License-Identifier: GPL-3.0-or-later

set(header_path "${BSMPT_SOURCE_DIR}/include/BSMPT/baryo_fhck")
set(header
    ${header_path}/TransportModel.h ${header_path}/BenchmarkModel.h
    ${header_path}/TransportEquations.h ${header_path}/LocalCalcEta.h
    ${header_path}/Kernels.h)

set(src TransportModel.cpp BenchmarkModel.cpp TransportEquations.cpp
        LocalCalcEta.cpp Kernels.cpp)

add_library(BaryoFHCK ${header} ${src})
target_link_libraries(
  BaryoFHCK
  PUBLIC Eigen3::Eigen
         GSL::gsl
         Minimizer
         Utility
         BounceSolution
         TransitionTracer
         GW
         VacuumProfile)

target_compile_definitions(
  BaryoFHCK
  PRIVATE KERNEL_DIR="${BSMPT_SOURCE_DIR}/include/BSMPT/baryo_fhck/kernels")
target_include_directories(BaryoFHCK PUBLIC ${BSMPT_SOURCE_DIR}/include)
target_compile_features(BaryoFHCK PUBLIC cxx_std_17)
