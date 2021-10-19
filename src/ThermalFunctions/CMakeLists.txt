# SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas Müller
#
# SPDX-License-Identifier: GPL-3.0-or-later

set(header_path "${BSMPT_SOURCE_DIR}/include/BSMPT/ThermalFunctions")
set(header
	${header_path}/ThermalFunctions.h
	${header_path}/NegativeBosonSpline.h
	${header_path}/thermalcoefficientcalculator.h
)
set(src
	ThermalFunctions.cpp
	thermalcoefficientcalculator.cpp
)

add_library(ThermalFunctions ${header} ${src})

target_link_libraries(ThermalFunctions PUBLIC GSL::gsl)
target_include_directories(ThermalFunctions PUBLIC ${BSMPT_SOURCE_DIR}/include)
target_compile_features(ThermalFunctions PUBLIC cxx_std_14)

