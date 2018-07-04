add_library(Minimizer
  Minfunc_gen.cpp
  MinimizeGSL.cpp
  Minimizer.cpp
)
target_link_libraries(Minimizer PUBLIC CMAES Eigen3::Eigen GSL::gsl)
target_compile_options(Minimizer PUBLIC ${OpenMP_CXX_FLAGS})
