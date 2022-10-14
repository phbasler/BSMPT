list(APPEND CMAKE_MODULE_PATH ${CMAKE_BINARY_DIR})
include(conan)
set(ConanPkgList eigen/3.4.0 gsl/2.7 boost/1.78.0 nlohmann_json/3.10.5)

if(BSMPT_IS_TOPLEVEL AND BUILD_TESTING)
  set(ConanPkgList ${ConanPkgList} catch2/3.1.0)
endif()

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  set(ConanPkgList ${ConanPkgList} benchmark/1.6.1)
endif()
if(UseNLopt)
  set(ConanPkgList ${ConanPkgList} nlopt/2.7.1)
endif(UseNLopt)
conan_cmake_configure(REQUIRES ${ConanPkgList} GENERATORS cmake_find_package)

conan_cmake_autodetect(settings)

conan_cmake_install(
  PATH_OR_REFERENCE
  .
  BUILD
  missing
  REMOTE
  conancenter
  SETTINGS
  ${settings})
