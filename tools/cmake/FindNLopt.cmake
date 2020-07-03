# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

#.rst:
# FindNLopt
# ----------
#
# Find NLopt include dirs and libraries.
#
# This module sets the following variables:
#
# ::
#
#   NLopt_FOUND - set to true if the library is found
#   NLopt_INCLUDE_DIRS - list of required include directories
#   NLopt_LIBRARIES - list of libraries to be linked
#

find_package(PkgConfig)
pkg_check_modules(PC_NLopt QUIET NLopt)

find_path(NLopt_INCLUDE_DIR
  NAMES nlopt.hpp
  PATHS ${PC_NLopt_INCLUDE_DIRS}
)

find_library(NLopt_LIBRARY
  NAMES nlopt nlopt_cxx
  PATHS ${PC_NLopt_LIBRARY_DIRS}
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(NLopt
  FOUND_VAR NLopt_FOUND
  REQUIRED_VARS
  NLopt_LIBRARY
  NLopt_INCLUDE_DIR
  VERSION_VAR NLopt_VERSION
)

if(NLopt_FOUND)
  set(NLopt_LIBRARIES ${NLopt_LIBRARY})
  set(NLopt_INCLUDE_DIRS ${NLopt_INCLUDE_DIR})
  set(NLopt_DEFINITIONS ${PC_NLopt_CFLAGS_OTHER})
endif()

mark_as_advanced(
  NLopt_INCLUDE_DIR
  NLopt_LIBRARY
)
