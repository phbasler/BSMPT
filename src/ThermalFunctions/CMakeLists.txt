set(header_path "${CMAKE_SOURCE_DIR}/include/BSMPT/ThermalFunctions")
set(header
	${header_path}/ThermalFunctions.h
	${header_path}/NegativeBosonSpline.h
)
set(src
	ThermalFunctions.cpp
)

add_library(ThermalFunctions ${header} ${src})

target_link_libraries(ThermalFunctions PUBLIC GSL::gsl)
target_include_directories(ThermalFunctions PUBLIC ${BSMPT_SOURCE_DIR}/include)
target_compile_features(ThermalFunctions PUBLIC cxx_std_14)

