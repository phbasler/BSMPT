# SPDX-FileCopyrightText: 2024 Lisa Biermann, Margarete Mühlleitner, Rui Santos,
# João Viana
#
# SPDX-License-Identifier: GPL-3.0-or-later

set(suffix "include/BSMPT/minimum_tracer")
set(header_path "${BSMPT_SOURCE_DIR}/${suffix}")
set(header ${header_path}/minimum_tracer.h)

set(src minimum_tracer.cpp)

add_library(MinimumTracer STATIC ${header} ${src})
target_link_libraries(MinimumTracer PUBLIC Eigen3::Eigen GSL::gsl Minimizer
                                           Utility)

target_include_directories(MinimumTracer PUBLIC ${BSMPT_SOURCE_DIR}/include)
target_compile_features(MinimumTracer PUBLIC cxx_std_17)

set_property(TARGET MinimumTracer PROPERTY PUBLIC_HEADER ${header})

install(TARGETS MinimumTracer PUBLIC_HEADER DESTINATION ${suffix})
