# SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
# Müller
#
# SPDX-License-Identifier: GPL-3.0-or-later

set(header_path "${BSMPT_SOURCE_DIR}/include/BSMPT/utility")
set(header ${header_path}/utility.h ${header_path}/Logger.h
           ${header_path}/parser.h ${header_path}/const_velocity_spline.h ${header_path}/NumericalDerivatives.h)
set(src utility.cpp Logger.cpp parser.cpp const_velocity_spline.cpp NumericalDerivatives.cpp)
add_library(Utility ${header} ${src})
target_include_directories(Utility PUBLIC ${BSMPT_SOURCE_DIR}/include
                                          ${BSMPT_BINARY_DIR}/include)
if(Boost_FOUND)
  target_link_libraries(Utility PUBLIC Boost::boost ASCIIPlotter Spline)
endif()
target_compile_features(Utility PUBLIC cxx_std_17)

if(nlohmann_json_FOUND)
  target_link_libraries(Utility PRIVATE nlohmann_json::nlohmann_json
                                        ASCIIPlotter Spline)
endif()
target_compile_features(Utility PUBLIC cxx_std_14)
