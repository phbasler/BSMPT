# SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
# Müller
#
# SPDX-License-Identifier: GPL-3.0-or-later

cmake_minimum_required(VERSION 3.23)
set(BSMPT_VERSION 3.1.3)
project(
  BSMPT
  VERSION ${BSMPT_VERSION}
  LANGUAGES C CXX
  DESCRIPTION
    "BSMPT - Beyond the Standard Model Phase Transitions : A C++ package for the computation of the EWPT in BSM models"
)

# -------------------------------- options     --------------------------------

option(UseLibCMAES "Use libcmaes to minimize the potential" ON)
option(UseNLopt "Use NLopt to minimize the potential" ON)
option(MakeAdditionalTesting "Build the executables used for further testing"
       OFF)
option(BSMPTCompileBaryo "Compile Baryogenesis" OFF)
option(EnableCoverage "Enable code coverage" OFF)
option(BSMPTUseVectorization "Disable vectorization" ON)
option(BSMPTBuildExecutables "Build the executables" ON)

set(BSMPT_IS_TOPLEVEL NO)
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  set(BSMPT_IS_TOPLEVEL YES)
endif()

if(NOT BSMPT_IS_TOPLEVEL)
  set(EnableCoverage OFF)
endif()

if(EnableCoverage)
  # The compiler flags don't mix well, so we turn it off in this scenario
  set(BSMPTUseVectorization OFF)
endif()

if(BSMPT_IS_TOPLEVEL)
  include(CTest)
endif()

# -------------------------------- setup cmake --------------------------------
if(POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW) # use <packagename>_ROOT as search path
endif(POLICY CMP0074)

list(APPEND CMAKE_MODULE_PATH "${BSMPT_SOURCE_DIR}/tools/cmake")

include(CompilerOptions)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif(NOT CMAKE_BUILD_TYPE)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin")

include(LTO)

find_package(PkgConfig)
find_package(Eigen3 CONFIG REQUIRED)
find_package(Boost CONFIG)
find_package(GSL CONFIG REQUIRED)
find_package(nlohmann_json CONFIG)

find_package(Threads REQUIRED)

if(NOT Boost_FOUND)
  set(BSMPTCompileBaryo OFF)
endif()

set(CodeCoverageExcludesFromOtherPkgs "")
if(NOT "${EIGEN3_ROOT_DIR}" STREQUAL "")
  list(APPEND CodeCoverageExcludesFromOtherPkgs "${EIGEN3_ROOT_DIR}/*")
endif()

if(NOT "${GSL_ROOT_DIR}" STREQUAL "")
  list(APPEND CodeCoverageExcludesFromOtherPkgs "${GSL_ROOT_DIR}/*")
endif()

if(Boost_FOUND)
  if(NOT "${Boost_DIR}" STREQUAL "")
    list(APPEND CodeCoverageExcludesFromOtherPkgs "${Boost_DIR}/*")
  endif()
  if(NOT "${Boost_ROOT_DIR}" STREQUAL "")
    list(APPEND CodeCoverageExcludesFromOtherPkgs "${Boost_ROOT_DIR}/*")
  endif()
endif()

if(nlohmann_json_FOUND)
  if(NOT "${nlohmann_json}" STREQUAL "")
    list(APPEND CodeCoverageExcludesFromOtherPkgs "${nlohmann_json}/*")
  endif()
endif()

set(EXPORT_CMAES FALSE)
include(FetchContent)
if(UseLibCMAES)
  find_package(cmaes CONFIG REQUIRED)
endif(UseLibCMAES)

if(UseNLopt)
  find_package(NLopt CONFIG QUIET)
  if(NLopt_FOUND)
    if(NOT "${NLopt_ROOT_DIR}" STREQUAL "")
      set(CodeCoverageExcludesFromOtherPkgs ${CodeCoverageExcludesFromOtherPkgs}
                                            "${NLopt_ROOT_DIR}/*")
    endif()
    if(NOT "${NLopt_DIR}" STREQUAL "")
      set(CodeCoverageExcludesFromOtherPkgs ${CodeCoverageExcludesFromOtherPkgs}
                                            "${NLopt_DIR}/*")
    endif()
  endif()

endif(UseNLopt)

if(NOT cmaes_FOUND AND NOT NLopt_FOUND)
  message(FATAL_ERROR "You need at least libcmaes or NLopt.")
endif()

# -------------------- Coverage
# --------------------------------------------------

if(EnableCoverage)
  add_compile_options(-coverage)
  add_link_options(-coverage)
  include(CodeCoverage)
  append_coverage_compiler_flags()
  # set(LCOV_ARGS "--no-external")
  set(COVERAGE_EXCLUDES
      "/usr/include/*"
      "/usr/lib/*"
      "~/.conan2/*"
      "~/.conan/*"
      "${PROJECT_SOURCE_DIR}/tools/*"
      ${CodeCoverageExcludesFromOtherPkgs}
      "${PROJECT_SOURCE_DIR}/tests/*"
      "tests/*"
      "${PROJECT_SOURCE_DIR}/src/prog/*"
      "${PROJECT_SOURCE_DIR}/src/utility/*"
      "${PROJECT_SOURCE_DIR}/include/BSMPT/utility/Logger.h"
      "${PROJECT_SOURCE_DIR}/src/Kfactors/Kfactors.cpp"
      "${PROJECT_SOURCE_DIR}/src/models/ClassPotentialOrigin_deprecated.cpp"
      "${PROJECT_SOURCE_DIR}/src/models/ClassTemplate.cpp"
      "${PROJECT_SOURCE_DIR}/src/WallThickness/WallThickness1D.cpp"
      "${PROJECT_SOURCE_DIR}/src/standalone/*"
      "${PROJECT_SOURCE_DIR}/src/utility/spline/spline.cpp"
      "${PROJECT_SOURCE_DIR}/src/utility/asciiplotter/asciiplotter.cpp")
  setup_target_for_coverage_lcov(
    NAME
    coverage # New target name
    EXECUTABLE
    ctest
    -j
    ${n_cores} # Executable in PROJECT_BINARY_DIR
    --output-on-failure
    BASE_DIRECTORY
    "${PROJECT_SOURCE_DIR}"
    DEPENDS
    GenericTests)
endif(EnableCoverage)

find_lto(CXX)

add_subdirectory(src)
add_subdirectory(standalone)

# ------------------------ Doxygen --------------------------------------------
if(BSMPT_IS_TOPLEVEL)
  find_package(Doxygen)
  if(DOXYGEN_FOUND)

    set(DOXYGEN_USE_MATHJAX YES)
    set(DOXYGEN_STRIP_FROM_INC_PATH "${CMAKE_CURRENT_SOURCE_DIR}/include")
    set(DOXYGEN_EXCLUDE
        ${CMAKE_CURRENT_SOURCE_DIR}/src/Kfactors/Kfunctions_grid_Part1.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/Kfactors/Kfunctions_grid_Part2.cpp)
    set(DOXYGEN_EXCLUDE_PATTERNS ${CMAKE_CURRENT_SOURCE_DIR}/src/*/README.md)
    set(DOXYGEN_PROJECT_BRIEF ${CMAKE_PROJECT_DESCRIPTION})
    set(DOXYGEN_EXTRACT_PRIVATE YES)
    set(DOXYGEN_GENERATE_TREEVIEW YES)
    set(DOXYGEN_DISTRIBUTE_GROUP_DOC YES)
    set(DOXYGEN_WARN_IF_UNDOCUMENTED YES)
    set(DOXYGEN_WARN_IF_DOC_ERROR YES)
    set(DOXYGEN_USE_MDFILE_AS_MAINPAGE ${CMAKE_CURRENT_SOURCE_DIR}/README.md)

    doxygen_add_docs(
      doc
      "${CMAKE_CURRENT_SOURCE_DIR}/include/"
      "${CMAKE_CURRENT_SOURCE_DIR}/src/"
      "${CMAKE_CURRENT_SOURCE_DIR}/README.md"
      "${CMAKE_CURRENT_SOURCE_DIR}/Changelog.md"
      "${CMAKE_CURRENT_SOURCE_DIR}/tests/")

  else()
    message(
      "Doxygen need to be installed to generate the doxygen documentation")
  endif()

endif(BSMPT_IS_TOPLEVEL)

# -------------------- Testing
# ---------------------------------------------------

include(CTest)
if(BSMPT_IS_TOPLEVEL AND BUILD_TESTING)
  add_subdirectory(tests)
endif()
