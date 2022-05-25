find_package(OpenMP REQUIRED)
find_package(libcmaes 0.10 QUIET)
  if(NOT libcmaes_FOUND)
    set(BSMPT_CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
    set(BSMPT_CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}")
    set(BSMPT_CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")
    set(EXPORT_CMAES TRUE)
    FetchContent_Declare(
      libcmaes
      GIT_REPOSITORY https://github.com/CMA-ES/libcmaes.git
      GIT_TAG v0.10)
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

    ## WORKAROUND for https://github.com/phbasler/BSMPT/issues/72 until the PR is fixed
    string(REPLACE " " ";" REPLACED_FLAGS ${BSMPT_CMAKE_CXX_FLAGS})
    string(REPLACE " " ";" REPLACED_FLAGS_RELEASE ${BSMPT_CMAKE_CXX_FLAGS_RELEASE})
    string(REPLACE " " ";" REPLACED_FLAGS_DEBUG ${BSMPT_CMAKE_CXX_FLAGS_DEBUG})
    target_compile_options(cmaes PUBLIC ${REPLACED_FLAGS})
    if(CMAKE_BUILD_TYPE MATCHES RELEASE)
      target_compile_options(cmaes PUBLIC ${REPLACED_FLAGS_RELEASE})
    endif()
    if(CMAKE_BUILD_TYPE MATCHES DEBUG)
      target_compile_options(cmaes PUBLIC ${REPLACED_FLAGS_DEBUG})
    endif()
    set(CMAKE_CXX_FLAGS "${BSMPT_CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS}")
    set(CMAKE_CXX_FLAGS_DEBUG "${BSMPT_CMAKE_CXX_FLAGS_DEBUG} ${CMAKE_CXX_FLAGS_DEBUG}")
    set(CMAKE_CXX_FLAGS_RELEASE "${BSMPT_CMAKE_CXX_FLAGS_RELEASE} ${CMAKE_CXX_FLAGS_RELEASE}")
    ## WORKAROUND end

    set(libcmaes_FOUND TRUE)
  else()
    set(CodeCoverageExcludesFromOtherPkgs
      ${CodeCoverageExcludesFromOtherPkgs}
      "${libcmaes_ROOT_DIR}/*"
      )
  endif()
