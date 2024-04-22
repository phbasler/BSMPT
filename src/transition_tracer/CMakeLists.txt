# SPDX-FileCopyrightText: 2024 Lisa Biermann, Margarete Mühlleitner, Rui
# Santos, João Viana
#
# SPDX-License-Identifier: GPL-3.0-or-later

set(header_path "${BSMPT_SOURCE_DIR}/include/BSMPT/transition_tracer")
set(header
  ${header_path}/transition_tracer.h
)

set(src
  transition_tracer.cpp
)

add_library(TransitionTracer ${header} ${src})
target_link_libraries(TransitionTracer PUBLIC MinimumTracer BounceSolution GW)
target_include_directories(TransitionTracer PUBLIC ${BSMPT_SOURCE_DIR}/include)
target_compile_features(TransitionTracer PUBLIC cxx_std_14)
