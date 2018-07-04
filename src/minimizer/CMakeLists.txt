add_library(Minimizer
  Minfunc_gen.cpp
  MinimizeGSL.cpp
  Minimizer.cpp
)
target_include_directories( Minimizer PUBLIC ${EIGEN3_INCLUDE_DIR})
target_link_libraries(Minimizer PRIVATE CMAES)
