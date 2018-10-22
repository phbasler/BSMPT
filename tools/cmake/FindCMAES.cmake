if(DEFINED ENV{CMAES_ROOT})
  set(CMAES_ROOT $ENV{CMAES_ROOT})
endif(DEFINED ENV{CMAES_ROOT})

find_library( CMAES_LIBRARY
  NAMES cmaes
  HINTS ${CMAES_ROOT}/lib ${CMAES_ROOT}/lib64 ${CMAKE_SOURCE_DIR}/tools/lib ${CMAKE_SOURCE_DIR}/tools/lib64
  DOC "path to cmaes library"
)

find_path( CMAES_INCLUDE_DIR
  NAMES libcmaes/cmaes.h
  HINTS ${CMAES_ROOT}/include/ ${CMAKE_SOURCE_DIR}/tools/include
  DOC "cmaes include directory"
)

if(NOT CMAES_LIBRARY OR NOT CMAES_INCLUDE_DIR)
  MESSAGE(STATUS "libCMAES not found, downloading...")
  file(MAKE_DIRECTORY ${CMAES_ROOT})
  file(
			DOWNLOAD
				https://github.com/beniz/libcmaes/archive/0.9.5.tar.gz
				${CMAKE_BINARY_DIR}/libcmaes.tar.gz
      EXPECTED_MD5
        ab89fde799f1e938ffd74fa66ab77153
		)
		execute_process(
			COMMAND
				${CMAKE_COMMAND} -E tar xzf ${CMAKE_BINARY_DIR}/libcmaes.tar.gz
		)
    execute_process(
			COMMAND ./autogen.sh
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/libcmaes-0.9.5
		)
    file(WRITE ${CMAKE_BINARY_DIR}/libcmaes-0.9.5/cmaes_export.h "#define CMAES_EXPORT")
    execute_process(
      COMMAND ./configure
        CC=${CMAKE_C_COMPILER}
        CXX=${CMAKE_CXX_COMPILER}
        --with-eigen3-include=${EIGEN3_INCLUDE_DIR}
        --prefix=${CMAES_ROOT}
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/libcmaes-0.9.5
    )
    execute_process(
      COMMAND make
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/libcmaes-0.9.5
    )
    execute_process(
      COMMAND make install
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/libcmaes-0.9.5
    )
    find_library( CMAES_LIBRARY
      NAMES cmaes
      HINTS ${CMAES_ROOT}/lib ${CMAES_ROOT}/lib64
      DOC "path to cmaes library"
    )

    find_path( CMAES_INCLUDE_DIR
      NAMES libcmaes/cmaes.h
      HINTS ${CMAES_ROOT}/include/
      DOC "cmaes include directory"
    )
endif(NOT CMAES_LIBRARY OR NOT CMAES_INCLUDE_DIR)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args( CMAES
  FOUND_VAR
    CMAES_FOUND
  REQUIRED_VARS
    CMAES_LIBRARY
    CMAES_INCLUDE_DIR
)

if(CMAES_FOUND AND NOT TARGET CMAES)
  add_library(CMAES UNKNOWN IMPORTED)
  set_property(
    TARGET CMAES
    PROPERTY IMPORTED_LOCATION ${CMAES_LIBRARY}
  )
  set_property(
    TARGET CMAES
    PROPERTY INTERFACE_INCLUDE_DIRECTORIES
      ${CMAES_INCLUDE_DIR}
  )
endif(CMAES_FOUND AND NOT TARGET CMAES)
