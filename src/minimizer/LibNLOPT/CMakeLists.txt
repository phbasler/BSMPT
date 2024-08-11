# SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
# Müller
#
# SPDX-License-Identifier: GPL-3.0-or-later

set(suffix "include/BSMPT/minimizer/LibNLOPT")
set(header_path "${BSMPT_SOURCE_DIR}/${suffix}")
set(header ${header_path}/MinimizeNLOPT.h)

set(src MinimizeNLOPT.cpp)

add_library(Minimizer_NLOPT INTERFACE)
target_sources(Minimizer_NLOPT INTERFACE ${header} ${src})
target_link_libraries(Minimizer_NLOPT INTERFACE Minimizer)

if(TARGET NLopt::nlopt)
  target_link_libraries(Minimizer_NLOPT INTERFACE NLopt::nlopt)
else()
  target_link_libraries(Minimizer_NLOPT INTERFACE NLopt::nlopt_cxx)
endif()

target_include_directories(
  Minimizer_NLOPT INTERFACE ${BSMPT_SOURCE_DIR}/include
                            ${BSMPT_BINARY_DIR}/include)
target_include_directories(Minimizer_NLOPT
                           INTERFACE ${BSMPT_SOURCE_DIR}/include)
target_compile_features(Minimizer_NLOPT INTERFACE cxx_std_17)

set_property(TARGET Minimizer_NLOPT PROPERTY PUBLIC_HEADER ${header})

install(TARGETS Minimizer_NLOPT PUBLIC_HEADER DESTINATION ${suffix})
