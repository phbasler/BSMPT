set(header_path "${CMAKE_SOURCE_DIR}/include/BSMPT/WallThickness")
set(header
	${header_path}/WallThicknessLib.h
)
set(src
	WallThicknessLib.cpp
)

add_library(LibWallThickness ${header} ${src})
target_link_libraries(LibWallThickness PUBLIC Eigen3::Eigen GSL::gsl Minimizer Boost::boost Threads::Threads)
target_include_directories(LibWallThickness PUBLIC ${BSMPT_SOURCE_DIR}/include)
target_compile_features(LibWallThickness PUBLIC cxx_std_14)

