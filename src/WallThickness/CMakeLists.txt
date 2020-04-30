set(header_path "${CMAKE_SOURCE_DIR}/include/BSMPT/WallThickness")
set(header
	${header_path}/WallThicknessLib.h
)
set(src
	WallThicknessLib.cpp
)

add_library(lib_WallThickness ${header} ${src})
target_link_libraries(lib_WallThickness PUBLIC Eigen3::Eigen GSL::gsl Minimizer Boost::boost)
target_include_directories(lib_WallThickness PUBLIC ${BSMPT_SOURCE_DIR}/include)
target_compile_features(lib_WallThickness PUBLIC cxx_std_14)

