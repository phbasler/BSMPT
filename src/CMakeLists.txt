# SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
# Müller
#
# SPDX-License-Identifier: GPL-3.0-or-later

configure_file(${BSMPT_SOURCE_DIR}/include/BSMPT/config.h.in
               ${BSMPT_BINARY_DIR}/include/BSMPT/config.h)

add_library(BSMPT_Config INTERFACE)
set(header_path ${BSMPT_BINARY_DIR}/include/BSMPT)
target_include_directories(
  BSMPT_Config
  INTERFACE $<BUILD_INTERFACE:${BSMPT_SOURCE_DIR}/include>
            $<BUILD_INTERFACE:${BSMPT_BINARY_DIR}/include>
            $<INSTALL_INTERFACE:include>)
set(header ${header_path}/config.h)
set_property(TARGET BSMPT_Config PROPERTY PUBLIC_HEADER ${header})

install(TARGETS BSMPT_Config PUBLIC_HEADER DESTINATION include/BSMPT)

add_subdirectory(minimizer)
add_subdirectory(models)
add_subdirectory(models/modeltests)
if(BSMPTCompileBaryo)
  add_subdirectory(baryo_calculation)
  add_subdirectory(WallThickness)
  add_subdirectory(Kfactors)
endif(BSMPTCompileBaryo)
add_subdirectory(ThermalFunctions)

add_subdirectory(minimum_tracer)
add_subdirectory(transition_tracer)
add_subdirectory(bounce_solution)
add_subdirectory(gravitational_waves)

# ==============================Utility==================================

add_subdirectory(utility)
add_subdirectory(utility/asciiplotter)
add_subdirectory(utility/spline)

# ==============================EXE======================================

if(BSMPTBuildExecutables)

  add_executable(BSMPT prog/BSMPT.cpp)
  target_link_libraries(BSMPT Minimizer Models Utility)
  target_compile_features(BSMPT PUBLIC cxx_std_17)

  add_executable(VEVEVO prog/VEVEVO.cpp)
  target_link_libraries(VEVEVO Minimizer Models Utility)
  target_compile_features(VEVEVO PUBLIC cxx_std_17)

  add_executable(NLOVEV prog/NLOVEV.cpp)
  target_link_libraries(NLOVEV Minimizer Models Utility)
  target_compile_features(NLOVEV PUBLIC cxx_std_17)

  add_executable(CalcCT prog/CalcCT.cpp)
  target_link_libraries(CalcCT Models Utility)
  target_compile_features(CalcCT PUBLIC cxx_std_17)

  add_executable(TripleHiggsCouplingsNLO prog/TripleHiggsNLO.cpp)
  target_link_libraries(TripleHiggsCouplingsNLO Models Utility)
  target_compile_features(TripleHiggsCouplingsNLO PUBLIC cxx_std_17)

  add_executable(Test prog/Test.cpp)
  target_link_libraries(Test Minimizer Models Utility ModelTests)
  target_compile_features(Test PUBLIC cxx_std_17)

  add_executable(CalcGW prog/CalcGW.cpp)
  target_link_libraries(CalcGW Utility GW MinimumTracer TransitionTracer)
  target_compile_features(CalcGW PUBLIC cxx_std_17)

  add_executable(CalcTemps prog/CalcTemps.cpp)
  target_link_libraries(CalcTemps Utility GW MinimumTracer TransitionTracer)
  target_compile_features(CalcTemps PUBLIC cxx_std_17)

  add_executable(MinimaTracer prog/MinimaTracer.cpp)
  target_link_libraries(MinimaTracer Minimizer Models Utility MinimumTracer
                        TransitionTracer)
  target_compile_features(MinimaTracer PUBLIC cxx_std_17)

  add_executable(PotPlotter prog/PotPlotter.cpp)
  target_link_libraries(PotPlotter Models Utility Minimizer MinimumTracer)
  target_compile_features(PotPlotter PUBLIC cxx_std_17)

  if(BSMPTCompileBaryo)
    add_executable(CalculateEWBG prog/CalculateEWBG.cpp)
    target_link_libraries(CalculateEWBG Minimizer Models Baryo Utility)
    target_compile_features(CalculateEWBG PUBLIC cxx_std_17)
  else(BSMPTCompileBaryo)
    message("CalculateEWBG will not compile as Boost is missing")
  endif(BSMPTCompileBaryo)

  if(BSMPTCompileBaryo)
    add_executable(PlotEWBG_vw prog/PlotEWBG_vw.cpp)
    target_link_libraries(PlotEWBG_vw Minimizer Models Baryo Utility)
    target_compile_features(PlotEWBG_vw PUBLIC cxx_std_17)
  else(BSMPTCompileBaryo)
    message("PlotEWBG_vw will not compile as Boost is missing")
  endif(BSMPTCompileBaryo)

  if(BSMPTCompileBaryo)
    add_executable(PlotEWBG_nL prog/PlotEWBG_nL.cpp)
    target_link_libraries(PlotEWBG_nL Minimizer Models Baryo Utility)
    target_compile_features(PlotEWBG_nL PUBLIC cxx_std_17)
  else(BSMPTCompileBaryo)
    message("PlotEWBG_nL will not compile as Boost is missing")
  endif(BSMPTCompileBaryo)

  # ======================= Additional testing
  # ====================================

  if(MakeAdditionalTesting)

    if(BSMPTCompileBaryo)
      add_executable(CreateMuGrid prog/CreateMuGrid.cpp)
      target_link_libraries(CreateMuGrid Minimizer Models Baryo Utility)
      target_compile_features(CreateMuGrid PUBLIC cxx_std_17)
    else(BSMPTCompileBaryo)
      message("CreateMuGrid will not compile as Boost is missing")
    endif(BSMPTCompileBaryo)

    if(BSMPTCompileBaryo)
      add_executable(Kfactor_grid prog/Kfactor_grid.cpp)
      target_link_libraries(Kfactor_grid Minimizer Models Kfactors Utility)
      target_compile_features(Kfactor_grid PUBLIC cxx_std_17)
    else(BSMPTCompileBaryo)
      message("Kfactor_grid will not compile as Boost is missing")
    endif(BSMPTCompileBaryo)

  endif(MakeAdditionalTesting)

endif(BSMPTBuildExecutables)
