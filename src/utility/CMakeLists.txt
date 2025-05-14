# SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
# Müller
#
# SPDX-License-Identifier: GPL-3.0-or-later

set(suffix "include/BSMPT/utility")
set(header_path "${BSMPT_SOURCE_DIR}/${suffix}")
set(header
    ${header_path}/utility.h
    ${header_path}/Logger.h
    ${header_path}/parser.h
    ${header_path}/const_velocity_spline.h
    ${header_path}/NumericalDerivatives.h
    ${header_path}/ModelIDs.h
    ${header_path}/settings.h)
set(src utility.cpp Logger.cpp parser.cpp const_velocity_spline.cpp
        NumericalDerivatives.cpp ModelIDs.cpp)
add_library(Utility ${header} ${src})
target_include_directories(Utility PUBLIC ${BSMPT_SOURCE_DIR}/include
                                          ${BSMPT_BINARY_DIR}/include)
if(Boost_FOUND)
  target_link_libraries(Utility PUBLIC Boost::boost)
endif()
target_compile_features(Utility PUBLIC cxx_std_17)

if(nlohmann_json_FOUND)
  target_link_libraries(Utility PRIVATE nlohmann_json::nlohmann_json)
endif()

target_link_libraries(Utility PUBLIC ASCIIPlotter Spline GSL::gsl BSMPT_Config)

set_property(TARGET Utility PROPERTY PUBLIC_HEADER ${header})

install(TARGETS Utility PUBLIC_HEADER DESTINATION ${suffix})
