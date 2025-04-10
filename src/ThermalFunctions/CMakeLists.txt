# SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
# Müller
#
# SPDX-License-Identifier: GPL-3.0-or-later

set(suffix "include/BSMPT/ThermalFunctions")
set(header_path "${BSMPT_SOURCE_DIR}/${suffix}")

set(header
    ${header_path}/ThermalFunctions.h ${header_path}/NegativeBosonSpline.h
    ${header_path}/thermalcoefficientcalculator.h)
set(src thermalcoefficientcalculator.cpp ThermalFunctions.cpp
        NegativeBosonSpline.cpp)

add_library(ThermalFunctions ${header} ${src})

target_link_libraries(ThermalFunctions PUBLIC GSL::gsl Utility)
target_include_directories(ThermalFunctions PUBLIC ${BSMPT_SOURCE_DIR}/include)
target_compile_features(ThermalFunctions PUBLIC cxx_std_17)

set_property(TARGET ThermalFunctions PROPERTY PUBLIC_HEADER ${header})

install(TARGETS ThermalFunctions PUBLIC_HEADER DESTINATION ${suffix})
