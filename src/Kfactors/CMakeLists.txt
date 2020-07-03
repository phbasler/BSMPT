set(header_path "${CMAKE_SOURCE_DIR}/include/BSMPT/Kfactors")
set(header
	${header_path}/Kfactors.h
	${header_path}/KfactorsinterpolatedGSL.h
	${header_path}/Kfactors_grid/KtildeInterpolation.h
	${header_path}/Kfactors_grid/Kfunctions_grid.h
)

set(src
	Kfunctions_grid_Part1.cpp
	Kfunctions_grid_Part2.cpp
	Kfactors.cpp
	KfactorsinterpolatedGSL.cpp
)

add_library(Kfactors ${header} ${src})

target_link_libraries(Kfactors PUBLIC Eigen3::Eigen GSL::gsl Boost::boost Minimizer)
target_compile_options(Kfactors PUBLIC ${OpenMP_CXX_FLAGS})
target_include_directories(Kfactors PUBLIC ${BSMPT_SOURCE_DIR}/include)
target_compile_features(Kfactors PUBLIC cxx_std_14)
