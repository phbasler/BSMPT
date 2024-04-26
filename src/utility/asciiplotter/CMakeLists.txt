# SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
# Müller
#
# SPDX-License-Identifier: GPL-3.0-or-later

set(header_path "${BSMPT_SOURCE_DIR}/include/BSMPT/utility/asciiplotter")
set(header ${header_path}/asciiplotter.h)
set(src asciiplotter.cpp)
add_library(ASCIIPlotter ${header} ${src})
target_include_directories(ASCIIPlotter PUBLIC ${BSMPT_SOURCE_DIR}/include
                                               ${BSMPT_BINARY_DIR}/include)
if(Boost_FOUND)
  target_link_libraries(ASCIIPlotter PUBLIC Boost::boost)
endif()
target_compile_features(ASCIIPlotter PUBLIC cxx_std_17)

if(nlohmann_json_FOUND)
  target_link_libraries(ASCIIPlotter PRIVATE nlohmann_json::nlohmann_json)
endif()
target_compile_features(ASCIIPlotter PUBLIC cxx_std_14)
