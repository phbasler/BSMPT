# SPDX-FileCopyrightText: 2024 Lisa Biermann, Margarete Mühlleitner, Rui Santos,
# João Viana
#
# SPDX-License-Identifier: GPL-3.0-or-later

set(suffix "include/BSMPT/gravitational_waves")
set(header_path "${BSMPT_SOURCE_DIR}/${suffix}")
set(header ${header_path}/gw.h)

set(src gw.cpp)

add_library(GW STATIC ${header} ${src})
target_link_libraries(GW PUBLIC Eigen3::Eigen GSL::gsl Minimizer Utility
                                BounceSolution)
target_include_directories(GW PUBLIC ${BSMPT_SOURCE_DIR}/include)
target_compile_features(GW PUBLIC cxx_std_17)

set_property(TARGET GW PROPERTY PUBLIC_HEADER ${header})

install(TARGETS GW PUBLIC_HEADER DESTINATION ${suffix})
