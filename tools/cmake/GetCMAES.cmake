if(CMAES_ROOT OR DEFINED ENV{CMAES_ROOT})
  find_package(CMAES REQUIRED)
else()
  include(FetchContent)
  FetchContent_Declare(
    libcmaes
    GIT_REPOSITORY https://github.com/beniz/libcmaes.git
    GIT_TAG 0.9.5)
  FetchContent_GetProperties(libcmaes)
  if(NOT libcmaes_POPULATED)
    FetchContent_Populate(libcmaes)
    message(STATUS "building libcmaes...")
    execute_process(COMMAND ./autogen.sh
                    WORKING_DIRECTORY ${libcmaes_SOURCE_DIR})
    file(WRITE ${libcmaes_SOURCE_DIR}/cmaes_export.h "#define CMAES_EXPORT")
    execute_process(
      COMMAND ./configure CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER}
              --with-eigen3-include=${EIGEN3_INCLUDE_DIR}
              --prefix=${libcmaes_BINARY_DIR}
      WORKING_DIRECTORY ${libcmaes_SOURCE_DIR})
    execute_process(COMMAND make WORKING_DIRECTORY ${libcmaes_SOURCE_DIR})
    execute_process(COMMAND make install
                    WORKING_DIRECTORY ${libcmaes_SOURCE_DIR})
    set(CMAES_ROOT ${libcmaes_BINARY_DIR})
    find_package(CMAES REQUIRED)
    set(CMAES_DEPENDENCY "set(CMAES_ROOT ${CMAES_ROOT})")
  endif()
endif()
