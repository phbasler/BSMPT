add_library(Models
  ClassPotentialC2HDM.cpp
  ClassPotentialOrigin.cpp
  ClassPotentialR2HDM.cpp
  ClassPotentialRN2HDM.cpp
  ClassTemplate.cpp
  IncludeAllModels.cpp
)
target_include_directories( Models PUBLIC ${EIGEN3_INCLUDE_DIR})
target_link_libraries( Models PUBLIC CMAES GSL::gsl)
