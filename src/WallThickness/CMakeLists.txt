# SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas Müller
#
# SPDX-License-Identifier: GPL-3.0-or-later

set(header_path "${BSMPT_SOURCE_DIR}/include/BSMPT/WallThickness")
set(header
	${header_path}/WallThicknessLib.h
)
set(src
	WallThicknessLib.cpp
)

add_library(LibWallThickness ${header} ${src})
target_link_libraries(LibWallThickness PUBLIC Eigen3::Eigen GSL::gsl Minimizer Utility Boost::boost Threads::Threads)
target_include_directories(LibWallThickness PUBLIC ${BSMPT_SOURCE_DIR}/include)
target_compile_features(LibWallThickness PUBLIC cxx_std_14)

