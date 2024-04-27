set(header_path "${BSMPT_SOURCE_DIR}/include/BSMPT/utility")
set(header ${header_path}/utility.h ${header_path}/Logger.h
           ${header_path}/parser.h)
set(src utility.cpp Logger.cpp parser.cpp)
add_library(Utility ${header} ${src})
target_include_directories(Utility PUBLIC ${BSMPT_SOURCE_DIR}/include
                                          ${BSMPT_BINARY_DIR}/include)
if(Boost_FOUND)
  target_link_libraries(Utility PUBLIC Boost::boost)
endif()
target_compile_features(Utility PUBLIC cxx_std_17)

if(nlohmann_json_FOUND)
  target_link_libraries(Utility PRIVATE nlohmann_json::nlohmann_json)
endif()
target_compile_features(Utility PUBLIC cxx_std_17)
