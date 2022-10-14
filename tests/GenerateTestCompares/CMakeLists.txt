# SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
# Müller
#
# SPDX-License-Identifier: GPL-3.0-or-later

set(header C2HDM.h R2HDM.h RN2HDM.h CPINTHEDARK.h CXSM.h)

set(src C2HDM.cpp R2HDM.cpp RN2HDM.cpp CPINTHEDARK.cpp CXSM.cpp)

add_library(TestCompares ${header} ${src})
target_compile_features(TestCompares PUBLIC cxx_std_14)
target_link_libraries(TestCompares PRIVATE Minimizer)
target_include_directories(TestCompares PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
