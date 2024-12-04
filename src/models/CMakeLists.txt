# SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
# Müller
#
# SPDX-License-Identifier: GPL-3.0-or-later

set(header_path "${BSMPT_SOURCE_DIR}/include/BSMPT/models")
set(header
    ${header_path}/SMparam.h
    ${header_path}/ModelTestfunctions.h
    ${header_path}/IncludeAllModels.h
    ${header_path}/ClassPotentialOrigin.h
    ${header_path}/ClassPotentialC2HDM.h
    ${header_path}/ClassPotentialR2HDM.h
    ${header_path}/ClassPotentialN2HDM.h
    ${header_path}/ClassPotentialCxSM.h
    ${header_path}/ClassPotentialSM.h
    ${header_path}/ClassPotentialCPintheDark.h
    # Here you can add your new model
    ${header_path}/ClassTemplate.h)

set(src
    SMParam.cpp
    ModelTestfunctions.cpp
    IncludeAllModels.cpp
    ClassPotentialOrigin.cpp
    ClassPotentialOrigin_deprecated.cpp
    ClassPotentialC2HDM.cpp
    ClassPotentialR2HDM.cpp
    ClassPotentialN2HDM.cpp
    ClassPotentialCxSM.cpp
    ClassPotentialSM.cpp
    ClassPotentialCPintheDark.cpp
    # Here you can add your new model
    ClassTemplate.cpp)

add_library(Models ${header} ${src})
target_link_libraries(Models PUBLIC GSL::gsl Eigen3::Eigen
                                    ThermalFunctions Utility)
target_include_directories(Models PUBLIC ${BSMPT_SOURCE_DIR}/include)
target_compile_features(Models PUBLIC cxx_std_17)

# Include code-coverage settings: target_link_libraries(Models PUBLIC
# coverage_config)
