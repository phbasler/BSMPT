# SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas Müller
#
# SPDX-License-Identifier: GPL-3.0-or-later

cmake_minimum_required(VERSION 3.13)
project(
  BSMPT
  VERSION 2.3.3
  LANGUAGES C CXX
  DESCRIPTION "BSMPT - Beyond the Standard Model Phase Transitions : A C++ package for the computation of the EWPT in BSM models")

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

option(UseLibCMAES "Use libcmaes to minimize the potential" ON)
option(UseNLopt "Use NLopt to minimize the potential" ON)
option(MakeAdditionalTesting "Build the executables used for further testing" OFF)
option(BSMPTCompileBaryo "Compile Baryogenesis" ON)

include(LTO)

option(UseConan "Use conan to get dependencies" OFF)

if(UseConan)
	list(APPEND CMAKE_MODULE_PATH ${CMAKE_BINARY_DIR})
	include(conan)
	set(ConanPkgList eigen/3.4.0 
				     gsl/2.7				     
				     boost/1.78.0)

	if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
	        set(ConanPkgList ${ConanPkgList} benchmark/1.6.1)
	endif()
	if(UseNLopt)
		set(ConanPkgList ${ConanPkgList} nlopt/2.7.1)
	endif(UseNLopt)
	conan_cmake_configure(REQUIRES ${ConanPkgList}
                      GENERATORS cmake_find_package)
					  
	
	conan_cmake_autodetect(settings)

	conan_cmake_install(PATH_OR_REFERENCE .
                    BUILD missing
                    REMOTE conancenter
                    SETTINGS ${settings})
	
endif(UseConan)


find_package(PkgConfig)
find_package(Eigen3 3.3.0 REQUIRED)
find_package(Boost 1.66)
find_package(GSL 2.1 REQUIRED)



find_package(Threads REQUIRED)

if(NOT BOOST_FOUND)
  set(BSMPTCompileBaryo OFF)
endif()


set(CodeCoverageExcludesFromOtherPkgs
  "${EIGEN3_ROOT_DIR}/*"
  "${GSL_ROOT_DIR}/*")

if(Boost_FOUND)
  if(NOT "${Boost_DIR}" STREQUAL "")
    set(CodeCoverageExcludesFromOtherPkgs
      ${CodeCoverageExcludesFromOtherPkgs}
      "${Boost_DIR}/*"
      )
  endif()
  if(NOT "${Boost_ROOT_DIR}" STREQUAL "")
    set(CodeCoverageExcludesFromOtherPkgs
      ${CodeCoverageExcludesFromOtherPkgs}
      "${Boost_ROOT_DIR}/*"
      )
  endif()
endif()

set(EXPORT_CMAES FALSE)
include(FetchContent)
if(UseLibCMAES AND NOT WIN32)
	include(CMAES)
endif(UseLibCMAES AND NOT WIN32)

if(UseNLopt)
	find_package(NLopt QUIET)
  if(NLopt_FOUND)
    if(NOT "${NLopt_ROOT_DIR}" STREQUAL "")
      set(CodeCoverageExcludesFromOtherPkgs
        ${CodeCoverageExcludesFromOtherPkgs}
        "${NLopt_ROOT_DIR}/*"
        )
    endif()
    if(NOT "${NLopt_DIR}" STREQUAL "")
      set(CodeCoverageExcludesFromOtherPkgs
        ${CodeCoverageExcludesFromOtherPkgs}
        "${NLopt_DIR}/*"
        )
    endif()
  endif()

endif(UseNLopt)

if(NOT libcmaes_FOUND AND NOT NLopt_FOUND)
  message( FATAL_ERROR "You need at least libcmaes or NLopt." )
endif()

configure_file(include/BSMPT/config.h.in include/BSMPT/config.h)

option(EnableCoverage "Enable code coverage" OFF)
if(EnableCoverage)
  include(CodeCoverage)
  append_coverage_compiler_flags()
#  set(LCOV_ARGS "--no-external")
  set(COVERAGE_EXCLUDES
            "/usr/include/*"
            "/usr/lib/*"
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
            )
  SETUP_TARGET_FOR_COVERAGE_LCOV(
     NAME coverage                 # New target name
     EXECUTABLE ctest  -j ${n_cores} # Executable in PROJECT_BINARY_DIR
     BASE_DIRECTORY "${PROJECT_SOURCE_DIR}"
#     EXCLUDE "_deps/*" "${PROJECT_SOURCE_DIR}/tests/*" "${PROJECT_SOURCE_DIR}/src/prog/*"
   )
endif(EnableCoverage)


find_lto(CXX)


add_subdirectory(src)


# -------------------------- export as cmake package --------------------------
if(PROJECT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  include(CMakePackageConfigHelpers)
  write_basic_package_version_file(
    BSMPTConfigVersion.cmake
    VERSION ${PACKAGE_VERSION}
    COMPATIBILITY AnyNewerVersion)

  set(ExportList
    Models
    Minimizer
    ThermalFunctions
    Utility
    )


  if(BSMPTCompileBaryo)
    set(ExportList
      ${ExportList}
      Baryo
      LibWallThickness
      Kfactors
      )
  endif(BSMPTCompileBaryo)


  export(
    TARGETS ${ExportList}
    NAMESPACE BSMPT::
    FILE BSMPTTargets.cmake
    )

  if(EXPORT_CMAES)
    export(
      TARGETS cmaes
      NAMESPACE libcmaes::
      APPEND FILE BSMPTTargets.cmake
      )
  endif(EXPORT_CMAES)




  configure_file(tools/cmake/BSMPTConfig.cmake.in BSMPTConfig.cmake @ONLY)


  export(PACKAGE BSMPT)
endif()



# ------------------------ Doxygen --------------------------------------------

find_package(Doxygen)
if(DOXYGEN_FOUND AND CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)

  set(DOXYGEN_USE_MATHJAX YES)
  set(DOXYGEN_STRIP_FROM_INC_PATH "${CMAKE_CURRENT_SOURCE_DIR}/include")
  set(DOXYGEN_EXCLUDE
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Kfactors/Kfunctions_grid_Part1.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Kfactors/Kfunctions_grid_Part2.cpp
    )
  set(DOXYGEN_PROJECT_BRIEF ${CMAKE_PROJECT_DESCRIPTION})
  set(DOXYGEN_EXTRACT_PRIVATE YES)
  set(DOXYGEN_GENERATE_TREEVIEW YES)
  set(DOXYGEN_DISTRIBUTE_GROUP_DOC YES)
  set(DOXYGEN_WARN_IF_UNDOCUMENTED YES)
  set(DOXYGEN_WARN_IF_DOC_ERROR YES)
  set(DOXYGEN_USE_MDFILE_AS_MAINPAGE README.md)

  doxygen_add_docs(doc
                     "${CMAKE_CURRENT_SOURCE_DIR}/include/"
                     "${CMAKE_CURRENT_SOURCE_DIR}/src/"
                     "${CMAKE_CURRENT_SOURCE_DIR}/README.md"
                     "${CMAKE_CURRENT_SOURCE_DIR}/Changelog.md")

else ()
  message("Doxygen need to be installed to generate the doxygen documentation")
endif()



# -------------------- Testing ---------------------------------------------------
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    include(CTest)
endif()

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME AND BUILD_TESTING)
    add_subdirectory(tests)
endif()
