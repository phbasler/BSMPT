find_package(OpenMP REQUIRED)
find_package(libcmaes 0.10 QUIET)
  if(NOT libcmaes_FOUND)
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
    set(libcmaes_FOUND TRUE)
  else()
    set(CodeCoverageExcludesFromOtherPkgs
      ${CodeCoverageExcludesFromOtherPkgs}
      "${libcmaes_ROOT_DIR}/*"
      )
  endif()
