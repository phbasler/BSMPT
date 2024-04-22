# SPDX-FileCopyrightText: 2024 Lisa Biermann, Margarete Mühlleitner, Rui
# Santos, João Viana
#
# SPDX-License-Identifier: GPL-3.0-or-later

set(header_path "${BSMPT_SOURCE_DIR}/include/BSMPT/gravitational_waves")
set(header
  ${header_path}/gw.h
)

set(src
  gw.cpp
)

add_library(GW ${header} ${src})
target_link_libraries(GW PUBLIC Eigen3::Eigen GSL::gsl Minimizer Utility BounceSolution)
target_include_directories(GW PUBLIC ${BSMPT_SOURCE_DIR}/include)
target_compile_features(GW PUBLIC cxx_std_14)