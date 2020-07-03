# -----------------------------------------------------------
# Search for the libcmaes minimization library
#
# Takes hints from the CMAES_ROOT cmake variable and the CMAES_ROOT and
# CMAES_ROOT_DIR environment variables.
#
# If succesful this creates the CMAES::CMAES imported target.
# ------------------------------------------------------------

find_library(
  CMAES_LIBRARY
  NAMES cmaes
  DOC "path to cmaes library")
find_path(
  CMAES_INCLUDE_DIR
  NAMES libcmaes/cmaes.h
  DOC "cmaes include directory")

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CMAES FOUND_VAR CMAES_FOUND REQUIRED_VARS
                                  CMAES_LIBRARY CMAES_INCLUDE_DIR)

if(CMAES_FOUND AND NOT TARGET CMAES)
  add_library(CMAES::cmaes UNKNOWN IMPORTED)
  set_target_properties(CMAES::cmaes PROPERTIES IMPORTED_LOCATION
                                                ${CMAES_LIBRARY})
  set_target_properties(CMAES::cmaes PROPERTIES INTERFACE_INCLUDE_DIRECTORIES
                                                ${CMAES_INCLUDE_DIR})
endif(CMAES_FOUND AND NOT TARGET CMAES)

mark_as_advanced(CMAES_LIBRARY CMAES_INCLUDE_DIR)
