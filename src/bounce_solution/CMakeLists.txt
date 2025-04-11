# SPDX-FileCopyrightText: 2024 Lisa Biermann, Margarete Mühlleitner, Rui Santos,
# João Viana
#
# SPDX-License-Identifier: GPL-3.0-or-later

set(suffix "include/BSMPT/bounce_solution")
set(header_path "${BSMPT_SOURCE_DIR}/${suffix}")
set(header ${header_path}/bounce_solution.h ${header_path}/action_calculation.h
           ${header_path}/gstar.h)

set(src bounce_solution.cpp action_calculation.cpp gstar.cpp)

add_library(BounceSolution ${header} ${src})
target_link_libraries(BounceSolution PUBLIC Eigen3::Eigen GSL::gsl Minimizer
                                            Utility MinimumTracer)
target_include_directories(BounceSolution PUBLIC ${BSMPT_SOURCE_DIR}/include)
target_compile_features(BounceSolution PUBLIC cxx_std_17)

set_property(TARGET BounceSolution PROPERTY PUBLIC_HEADER ${header})

install(TARGETS BounceSolution PUBLIC_HEADER DESTINATION ${suffix})
