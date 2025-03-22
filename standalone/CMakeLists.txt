# SPDX-FileCopyrightText: 2024 Lisa Biermann, Margarete Mühlleitner, Rui Santos,
# João Viana
#
# SPDX-License-Identifier: GPL-3.0-or-later

if(BSMPTBuildExecutables)

  set(STANDALONE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
  file(GLOB SOURCE_FILES "*.cpp" CONFIURE_DEPENDS)

  foreach(file_path ${SOURCE_FILES})
    get_filename_component(filename "${file_path}" NAME)
    string(REPLACE ".cpp" "" filename "${filename}")
    add_executable("${filename}" "${file_path}")
    target_link_libraries(
      "${filename}"
      PRIVATE Minimizer
              MinimumTracer
              Models
              Utility
              TestCompares
              ThermalFunctions
              BounceSolution
              TransitionTracer
              GW)
    set_target_properties(
      "${filename}" PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                               "${CMAKE_BINARY_DIR}/bin/standalone")

    target_compile_features("${filename}" PUBLIC cxx_std_17)
  endforeach()

endif(BSMPTBuildExecutables)
