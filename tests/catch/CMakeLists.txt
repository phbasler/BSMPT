add_library(Catch catch.cpp)
target_link_libraries(Catch PRIVATE Utility)
target_include_directories(Catch INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_features(Catch PUBLIC cxx_std_14)

target_compile_options(
  Catch
  PUBLIC
    $<$<CXX_COMPILER_ID:MSVC>:/DCATCH_CONFIG_ENABLE_BENCHMARKING>
    $<$<OR:$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>>:-DCATCH_CONFIG_ENABLE_BENCHMARKING>

    $<$<CXX_COMPILER_ID:MSVC>:/DCATCH_CONFIG_CONSOLE_WIDTH=180>
    $<$<OR:$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>>:-DCATCH_CONFIG_CONSOLE_WIDTH=180>
)
