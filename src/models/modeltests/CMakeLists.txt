# SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
# Müller
#
# SPDX-License-Identifier: GPL-3.0-or-later

set(header_path "${BSMPT_SOURCE_DIR}/include/BSMPT/models/modeltests")
set(header ${header_path}/ModelTestfunctions.h)

set(src ModelTestfunctions.cpp)

add_library(ModelTests ${header} ${src})
target_link_libraries(ModelTests PUBLIC Models Minimizer)
target_include_directories(ModelTests PUBLIC ${BSMPT_SOURCE_DIR}/include)
target_compile_features(ModelTests PUBLIC cxx_std_17)

# Include code-coverage settings: target_link_libraries(ModelTests PUBLIC
# coverage_config)
