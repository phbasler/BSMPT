find_package(Doxygen)
if(DOXYGEN_FOUND)

  set(DOXYGEN_USE_MATHJAX YES)
  set(DOXYGEN_STRIP_FROM_INC_PATH "${CMAKE_CURRENT_SOURCE_DIR}/include")
  set(DOXYGEN_EXCLUDE
      ${CMAKE_CURRENT_SOURCE_DIR}/src/Kfactors/Kfunctions_grid_Part1.cpp
      ${CMAKE_CURRENT_SOURCE_DIR}/src/Kfactors/Kfunctions_grid_Part2.cpp)
  set(DOXYGEN_EXCLUDE_PATTERNS
      ${CMAKE_CURRENT_SOURCE_DIR}/src/*/README.md
      ${CMAKE_CURRENT_SOURCE_DIR}/src/Kfactors/Kfunctions_grid_Part1.cpp
      ${CMAKE_CURRENT_SOURCE_DIR}/src/Kfactors/Kfunctions_grid_Part2.cpp
      ${CMAKE_CURRENT_SOURCE_DIR}/src/*
      ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmarks/benchmark-ewpt-c2hdm.cpp # THIS ONE
      ${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmarks/benchmark-ewbg-c2hdm.cpp # THIS ONE
      ${CMAKE_CURRENT_SOURCE_DIR}/tests/unittests) # THIS ONE

  set(DOXYGEN_PROJECT_BRIEF ${CMAKE_PROJECT_DESCRIPTION})
  set(DOXYGEN_EXTRACT_PRIVATE YES)
  set(DOXYGEN_GENERATE_TREEVIEW YES)
  set(DOXYGEN_DISTRIBUTE_GROUP_DOC YES)
  set(DOXYGEN_WARN_IF_UNDOCUMENTED YES)
  set(DOXYGEN_WARN_IF_DOC_ERROR YES)
  set(DOXYGEN_USE_MDFILE_AS_MAINPAGE ${CMAKE_CURRENT_SOURCE_DIR}/README.md)
  set(DOXYGEN_GENERATE_XML YES)
  set(DOXYGEN_EXTRACT_ALL YES)
  set(DOXYGEN_EXTRACT_STATIC YES)
  set(DOXYGEN_EXTRACT_LOCAL_CLASSES YES)
  set(DOXYGEN_ENABLE_PREPROCESSING YES)
  set(DOXYGEN_MACRO_EXPANSION YES)

  doxygen_add_docs(
    doc
    "${CMAKE_CURRENT_SOURCE_DIR}/include/"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/"
    "${CMAKE_CURRENT_SOURCE_DIR}/README.md"
    "${CMAKE_CURRENT_SOURCE_DIR}/Changelog.md"
    "${CMAKE_CURRENT_SOURCE_DIR}/tests/")

else()
  message("Doxygen need to be installed to generate the doxygen documentation")
endif()
