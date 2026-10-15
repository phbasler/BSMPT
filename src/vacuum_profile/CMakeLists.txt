# SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
# Müller
#
# SPDX-License-Identifier: GPL-3.0-or-later

set(header_path "${BSMPT_SOURCE_DIR}/include/BSMPT/vacuum_profile")
set(header ${header_path}/difeq_vacuum_profile.h
           ${header_path}/vacuum_profile.h)

set(src vacuum_profile.cpp)

add_library(VacuumProfile ${header} ${src})
target_link_libraries(VacuumProfile PUBLIC Eigen3::Eigen GSL::gsl Minimizer
                                           Utility)
# target_compile_options(BaryoFHCK PUBLIC ${OpenMP_CXX_FLAGS})
target_include_directories(VacuumProfile PUBLIC ${BSMPT_SOURCE_DIR}/include)
target_compile_features(VacuumProfile PUBLIC cxx_std_17)
