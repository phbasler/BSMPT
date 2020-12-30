add_library(Catch catch.cpp)
target_include_directories(Catch INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

set(TESTDIR ${CMAKE_CURRENT_SOURCE_DIR})

file(GLOB SOURCE_FILES "Test-rn2hdm.cpp" CONFIURE_DEPENDS)
add_executable(GenericTests ${SOURCE_FILES})
target_link_libraries(GenericTests Catch Minimizer Models Utility)
target_compile_features(BSMPT PUBLIC cxx_std_14)

include(Catch)
catch_discover_tests(GenericTests)
