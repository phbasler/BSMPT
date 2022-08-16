set(header_path "${BSMPT_SOURCE_DIR}/include/BSMPT/utility")
set(header
  ${header_path}/utility.h
  ${header_path}/Logger.h)
set(src
  utility.cpp
  Logger.cpp)
add_library(Utility ${header} ${src})
target_include_directories(Utility PUBLIC ${BSMPT_SOURCE_DIR}/include ${BSMPT_BINARY_DIR}/include)
if(Boost_FOUND)
  target_link_libraries(Utility PUBLIC Boost::boost)
endif()
target_compile_features(Utility PUBLIC cxx_std_14)

