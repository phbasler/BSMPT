# SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
# Müller
#
# SPDX-License-Identifier: GPL-3.0-or-later

set(header_path "${BSMPT_SOURCE_DIR}/include/BSMPT/baryo_calculation")
set(header
    ${header_path}/transport_equations.h
    ${header_path}/Fluid_Type/gen_calc.h
    ${header_path}/Fluid_Type/gen_func_fluid.h
    ${header_path}/Fluid_Type/top_source.h
    ${header_path}/Fluid_Type/bot_source.h
    ${header_path}/Fluid_Type/tau_source.h
    ${header_path}/CalculateEtaInterface.h)

set(src
    transport_equations.cpp
    Fluid_Type/gen_calc.cpp
    Fluid_Type/gen_func_fluid.cpp
    Fluid_Type/top_source.cpp
    Fluid_Type/bot_source.cpp
    Fluid_Type/tau_source.cpp
    CalculateEtaInterface.cpp)

add_library(Baryo ${header} ${src})
target_link_libraries(
  Baryo
  PUBLIC Eigen3::Eigen
         GSL::gsl
         Minimizer
         Kfactors
         LibWallThickness
         Utility
         Boost::boost)
target_include_directories(Baryo PUBLIC ${BSMPT_SOURCE_DIR}/include)
target_compile_features(Baryo PUBLIC cxx_std_17)
