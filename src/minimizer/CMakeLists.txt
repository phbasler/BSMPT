set(header_path "${CMAKE_SOURCE_DIR}/include/BSMPT/minimizer")
set(header
	${header_path}/Minimizer.h
	${header_path}/MinimizePlane.h
	${header_path}/MinimizeGSL.h
)

set(src
	MinimizeGSL.cpp
	Minimizer.cpp
	MinimizePlane.cpp
)


if(libcmaes_FOUND)
	set(header
		${header}
		${header_path}/LibCMAES/MinimizeLibCMAES.h
		)
	set(src
		${src}
		LibCMAES/MinimizeLibCMAES.cpp
		)
endif()

if(NLopt_FOUND)
	set(header
		${header}
		${header_path}/LibNLOPT/MinimizeNLOPT.h
		)
	set(src
		${src}
		LibNLOPT/MinimizeNLOPT.cpp
		)
endif()



add_library(Minimizer ${header} ${src})
target_link_libraries(Minimizer PUBLIC Eigen3::Eigen GSL::gsl)
if(libcmaes_FOUND)
  target_link_libraries(Minimizer PRIVATE libcmaes::cmaes OpenMP::OpenMP_CXX)
endif()

if(NLopt_FOUND)
  target_link_libraries(Minimizer PUBLIC NLopt::nlopt)
endif()

target_include_directories(Minimizer PUBLIC ${BSMPT_SOURCE_DIR}/include
	                                         ${BSMPT_BINARY_DIR}/include)
target_include_directories(Minimizer PUBLIC ${BSMPT_SOURCE_DIR}/include)
target_compile_features(Minimizer PUBLIC cxx_std_14)

