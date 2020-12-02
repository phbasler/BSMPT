set(header_path "${CMAKE_SOURCE_DIR}/include/BSMPT/models")
set(header

	${header_path}/SMparam.h

	${header_path}/IncludeAllModels.h
	${header_path}/ClassPotentialOrigin.h
	${header_path}/ClassPotentialC2HDM.h
	${header_path}/ClassPotentialR2HDM.h
  ${header_path}/ClassPotentialRN2HDM.h
  ${header_path}/ClassPotentialCxSM.h
  ${header_path}/ClassPotentialJonasSUSY.h

  # Here you can add your new model
  ${header_path}/ClassTemplate.h
)

set(src
  IncludeAllModels.cpp
  ClassPotentialOrigin.cpp
  ClassPotentialOrigin_deprecated.cpp

  ClassPotentialC2HDM.cpp
  ClassPotentialR2HDM.cpp
  ClassPotentialRN2HDM.cpp
  ClassPotentialCxSM.cpp
  ClassPotentialJonasSUSY.cpp


  # Here you can add your new model
  ClassTemplate.cpp
)

add_library(Models ${header} ${src})
target_link_libraries( Models PUBLIC GSL::gsl Eigen3::Eigen Minimizer ThermalFunctions Utility)
target_include_directories(Models PUBLIC ${BSMPT_SOURCE_DIR}/include)
target_compile_features(Models PUBLIC cxx_std_14)

