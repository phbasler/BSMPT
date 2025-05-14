# SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
# Müller
#
# SPDX-License-Identifier: GPL-3.0-or-later

set(suffix "include/BSMPT/utility/spline")
set(header_path "${BSMPT_SOURCE_DIR}/${suffix}")

set(header ${header_path}/spline.h)
set(src spline.cpp)
add_library(Spline ${header} ${src})
target_include_directories(Spline PUBLIC ${BSMPT_SOURCE_DIR}/include
                                         ${BSMPT_BINARY_DIR}/include)
if(Boost_FOUND)
  target_link_libraries(Spline PUBLIC Boost::boost)
endif()
target_compile_features(Spline PUBLIC cxx_std_17)

if(nlohmann_json_FOUND)
  target_link_libraries(Spline PRIVATE nlohmann_json::nlohmann_json)
endif()

target_compile_definitions(Spline PRIVATE _USE_MATH_DEFINES)

set_property(TARGET Spline PROPERTY PUBLIC_HEADER ${header})

install(TARGETS Spline PUBLIC_HEADER DESTINATION ${suffix})
