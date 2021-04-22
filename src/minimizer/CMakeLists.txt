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
target_link_libraries(Minimizer PUBLIC Eigen3::Eigen GSL::gsl Threads::Threads)
if(libcmaes_FOUND)
  if(TARGET libcmaes::cmaes)
    target_link_libraries(Minimizer PRIVATE libcmaes::cmaes OpenMP::OpenMP_CXX)
  elseif(TARGET BSMPT::cmaes)
    target_link_libraries(Minimizer PRIVATE BSMPT::cmaes OpenMP::OpenMP_CXX)
    endif()

endif()


if(NLopt_FOUND)
  if(TARGET NLopt::nlopt)
    target_link_libraries(Minimizer PUBLIC NLopt::nlopt)
  else()
    target_link_libraries(Minimizer PUBLIC NLopt::nlopt_cxx)
  endif()
endif()

target_include_directories(Minimizer PUBLIC ${BSMPT_SOURCE_DIR}/include
	                                         ${BSMPT_BINARY_DIR}/include)
target_include_directories(Minimizer PUBLIC ${BSMPT_SOURCE_DIR}/include)
target_compile_features(Minimizer PUBLIC cxx_std_14)

