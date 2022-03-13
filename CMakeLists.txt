# SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas Müller
#
# SPDX-License-Identifier: GPL-3.0-or-later

cmake_minimum_required(VERSION 3.13)
project(
  BSMPT
  VERSION 2.3.3
  LANGUAGES CXX)

# -------------------------------- setup cmake --------------------------------
if(POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW) # use <packagename>_ROOT as search path
endif(POLICY CMP0074)

list(APPEND CMAKE_MODULE_PATH "${BSMPT_SOURCE_DIR}/tools/cmake")

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif(NOT CMAKE_BUILD_TYPE)

set(CMAKE_CXX_FLAGS_DEBUG
    "${CMAKE_CXX_FLAGS_DEBUG} -Wall  -DCOMPILEDEBUG=true"
)
set(CMAKE_CXX_FLAGS_RELEASE
    "${CMAKE_CXX_FLAGS_RELEASE} -Wall ")
	
if(CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS_DEBUG
       "${CMAKE_CXX_FLAGS_DEBUG} -Wextra -Wmissing-declarations -Wmissing-include-dirs -Wshadow -pedantic")
  set(CMAKE_CXX_FLAGS_RELEASE
    "${CMAKE_CXX_FLAGS_RELEASE} -O3 -pedantic -Wextra -Wshadow")
endif(CMAKE_COMPILER_IS_GNUCXX)

if (MSVC)
     set(CMAKE_CXX_FLAGS_DEBUG
       "${CMAKE_CXX_FLAGS_DEBUG} /permissive- /bigobj /w44101")

     set(CMAKE_CXX_FLAGS_RELEASE
        "${CMAKE_CXX_FLAGS_RELEASE} /permissive- /bigobj /w44101 /Ox")

     set(MSVC_DISABLED_WARNINGS_LIST
      "C4061" # enumerator 'identifier' in switch of enum 'enumeration' is not
              # explicitly handled by a case label
              # Disable this because it flags even when there is a default.
      "C4100" # 'exarg' : unreferenced formal parameter
      "C4127" # conditional expression is constant
      "C4200" # nonstandard extension used : zero-sized array in
              # struct/union.
      "C4204" # nonstandard extension used: non-constant aggregate initializer
      "C4221" # nonstandard extension used : 'identifier' : cannot be
              # initialized using address of automatic variable
      "C4242" # 'function' : conversion from 'int' to 'uint8_t',
              # possible loss of data
      "C4244" # 'function' : conversion from 'int' to 'uint8_t',
              # possible loss of data
      "C4245" # 'initializing' : conversion from 'long' to
              # 'unsigned long', signed/unsigned mismatch
      "C4267" # conversion from 'size_t' to 'int', possible loss of data
      "C4371" # layout of class may have changed from a previous version of the
              # compiler due to better packing of member '...'
      "C4388" # signed/unsigned mismatch
      "C4296" # '>=' : expression is always true
      "C4350" # behavior change: 'std::_Wrap_alloc...'
      "C4365" # '=' : conversion from 'size_t' to 'int',
              # signed/unsigned mismatch
      "C4389" # '!=' : signed/unsigned mismatch
      "C4464" # relative include path contains '..'
      "C4510" # 'argument' : default constructor could not be generated
      "C4571"
      "C4512" # 'argument' : assignment operator could not be generated
      "C4514" # 'function': unreferenced inline function has been removed
      "C4548" # expression before comma has no effect; expected expression with
              # side-effect" caused by FD_* macros.
      "C4610" # struct 'argument' can never be instantiated - user defined
              # constructor required.
      "C4623" # default constructor was implicitly defined as deleted
      "C4625" # copy constructor could not be generated because a base class
              # copy constructor is inaccessible or deleted
      "C4626" # assignment operator could not be generated because a base class
              # assignment operator is inaccessible or deleted
	  "C4643"
      "C4668" # 'symbol' is not defined as a preprocessor macro, replacing with
              # '0' for 'directives'
              # Disable this because GTest uses it everywhere.
      "C4706" # assignment within conditional expression
      "C4710" # 'function': function not inlined
      "C4711" # function 'function' selected for inline expansion
      "C4800" # 'int' : forcing value to bool 'true' or 'false'
              # (performance warning)
      "C4820" # 'bytes' bytes padding added after construct 'member_name'
	  "C4868"
	  "C4996"
      "C5026" # move constructor was implicitly defined as deleted
      "C5027" # move assignment operator was implicitly defined as deleted
      "C5031"
	  "C5039"
      "C5045"
      )
      string(REPLACE "C" " -wd" MSVC_DISABLED_WARNINGS_STR
                            ${MSVC_DISABLED_WARNINGS_LIST})

      set(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} ${MSVC_DISABLED_WARNINGS_STR}")
 endif()


set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin")

option(UseLibCMAES "Use libcmaes to minimize the potential" ON)
option(UseNLopt "Use NLopt to minimize the potential" ON)
option(MakeAdditionalTesting "Build the executables used for further testing" OFF)

include(LTO)

option(UseConan "Use conan to get dependencies" OFF)

if(UseConan)
	list(APPEND CMAKE_MODULE_PATH ${CMAKE_BINARY_DIR})
	include(conan)
	set(ConanPkgList eigen/3.4.0 
				     gsl/2.7				     
				     boost/1.78.0)
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
  find_package(OpenMP REQUIRED)
  find_package(libcmaes 0.10 QUIET)
    if(NOT libcmaes_FOUND)
      set(EXPORT_CMAES TRUE)
      FetchContent_Declare(
        libcmaes
        GIT_REPOSITORY https://github.com/CMA-ES/libcmaes.git
        GIT_TAG master)
      FetchContent_GetProperties(libcmaes)
      if(NOT libcmaes_POPULATED)
        option(LIBCMAES_BUILD_TESTS OFF)
        option(LIBCMAES_BUILD_PYTHON OFF)
        option(LIBCMAES_BUILD_EXAMPLES OFF)
        FetchContent_Populate(libcmaes)
        add_subdirectory(${libcmaes_SOURCE_DIR} ${libcmaes_BINARY_DIR})
        include (GenerateExportHeader)
        generate_export_header (cmaes EXPORT_FILE_NAME ${libcmaes_SOURCE_DIR}/include/libcmaes/cmaes_export.h)
        message(STATUS ${libcmaes_SOURCE_DIR})
        set(CodeCoverageExcludesFromOtherPkgs
          ${CodeCoverageExcludesFromOtherPkgs}
          "${libcmaes_SOURCE_DIR}/*"
          )
      endif()
      set(libcmaes_FOUND TRUE)
    else()
      set(CodeCoverageExcludesFromOtherPkgs
        ${CodeCoverageExcludesFromOtherPkgs}
        "${libcmaes_ROOT_DIR}/*"
        )
    endif()
endif(UseLibCMAES AND NOT WIN32)

if(UseNLopt)
	find_package(NLopt REQUIRED)
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


  if(Boost_FOUND)
    set(ExportList
      ${ExportList}
      Baryo
      LibWallThickness
      Kfactors
      )
  endif(Boost_FOUND)


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



# -------------------- Testing ---------------------------------------------------
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    include(CTest)
endif()

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME AND BUILD_TESTING)
    add_subdirectory(tests)
endif()
