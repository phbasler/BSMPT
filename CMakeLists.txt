cmake_minimum_required(VERSION 3.13)
project(
  BSMPT
  VERSION 1.9.92
  LANGUAGES C CXX)

# -------------------------------- setup cmake --------------------------------
if(POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW) # use <packagename>_ROOT as search path
endif(POLICY CMP0074)

list(APPEND CMAKE_MODULE_PATH "${BSMPT_SOURCE_DIR}/tools/cmake")

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "RELEASE")
endif(NOT CMAKE_BUILD_TYPE)

set(CMAKE_CXX_FLAGS_DEBUG
    "${CMAKE_CXX_FLAGS_DEBUG} -pedantic -Wall -Wextra -Wmissing-declarations -Wmissing-include-dirs -Wshadow  -DCOMPILEDEBUG=true"
)
set(CMAKE_CXX_FLAGS_RELEASE
    "${CMAKE_CXX_FLAGS_RELEASE} -pedantic -Wall -Wextra -Wshadow")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin")

option(UseLibCMAES "Use libcmaes to minimize the potential" ON)
option(UseNLopt "Use NLopt to minimize the potential" ON)
option(MakeAdditionalTesting "Build the executables used for further testing" OFF)

find_package(Eigen3 3.3.0 REQUIRED)
find_package(Boost 1.66)
find_package(GSL 2.1 REQUIRED)
if(UseLibCMAES AND NOT WIN32)
  find_package(OpenMP REQUIRED)
  include(GetCMAES)
endif(UseLibCMAES AND NOT WIN32)

if(UseNLopt)
  find_package(NLopt QUIET)
endif(UseNLopt)

if(NOT CMAES_FOUND AND NOT NLopt_FOUND)
  message( FATAL_ERROR "You need at least libcmaes or NLopt." )
endif()

configure_file(include/BSMPT/config.h.in include/BSMPT/config.h)


add_subdirectory(src)



# -------------------------- export as cmake package --------------------------
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  BSMPTConfigVersion.cmake
  VERSION ${PACKAGE_VERSION}
  COMPATIBILITY AnyNewerVersion)

export(
  TARGETS Models Minimizer ThermalFunctions Utility
  NAMESPACE BSMPT::
  FILE BSMPTTargets.cmake
  )

configure_file(tools/cmake/BSMPTConfig.cmake.in BSMPTConfig.cmake @ONLY)

if(PROJECT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  export(PACKAGE BSMPT)
endif()



# ------------------------ Doxygen --------------------------------------------

find_package(Doxygen)
if(DOXYGEN_FOUND)

  set(DOXYGEN_USE_MATHJAX YES)
  set(DOXYGEN_STRIP_FROM_INC_PATH "${CMAKE_CURRENT_SOURCE_DIR}/include")
  set(DOXYGEN_EXCLUDE
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Kfactors/Kfunctions_grid_Part1.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Kfactors/Kfunctions_grid_Part2.cpp
    )
  set(DOXYGEN_PROJECT_BRIEF "BSMPT - Beyond the Standard Model Phase Transitions : A C++ package for the computation of the EWPT in BSM models")
  set(DOXYGEN_EXTRACT_PRIVATE YES)
  set(DOXYGEN_GENERATE_TREEVIEW YES)
  set(DOXYGEN_DISTRIBUTE_GROUP_DOC YES)
  set(DOXYGEN_WARN_IF_UNDOCUMENTED YES)
  set(DOXYGEN_WARN_IF_DOC_ERROR YES)

  doxygen_add_docs(doc
                     "${CMAKE_CURRENT_SOURCE_DIR}/include/"
                     "${CMAKE_CURRENT_SOURCE_DIR}/src/"
                     "${CMAKE_CURRENT_SOURCE_DIR}/README.md"
                     "${CMAKE_CURRENT_SOURCE_DIR}/Changelog.md")

else (DOXYGEN_FOUND)
  message("Doxygen need to be installed to generate the doxygen documentation")
endif(DOXYGEN_FOUND)
