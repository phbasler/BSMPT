add_library(Minimizer
  Minfunc_gen.cpp
  MinimizeGSL.cpp
  Minimizer.cpp
)
target_link_libraries(Minimizer PUBLIC CMAES Eigen3::Eigen GSL::gsl OpenMP::OpenMP_CXX)
