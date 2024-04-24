add_library(CatchEventListener catcheventlistener.cpp)
target_link_libraries(CatchEventListener PRIVATE Utility Catch2::Catch2WithMain)
target_include_directories(CatchEventListener
                           INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_features(CatchEventListener PUBLIC cxx_std_17)
