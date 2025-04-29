# SPDX-FileCopyrightText: 2024 Lisa Biermann, Margarete Mühlleitner, Rui Santos,
# João Viana
#
# SPDX-License-Identifier: GPL-3.0-or-later

set(suffix "include/BSMPT/transition_tracer")
set(header_path "${BSMPT_SOURCE_DIR}/${suffix}")

set(header ${header_path}/transition_tracer.h)

set(src transition_tracer.cpp)

add_library(TransitionTracer ${header} ${src})
target_link_libraries(TransitionTracer PUBLIC MinimumTracer BounceSolution GW)
target_include_directories(TransitionTracer PUBLIC ${BSMPT_SOURCE_DIR}/include)
target_compile_features(TransitionTracer PUBLIC cxx_std_17)

set_property(TARGET TransitionTracer PROPERTY PUBLIC_HEADER ${header})

install(TARGETS TransitionTracer PUBLIC_HEADER DESTINATION ${suffix})
