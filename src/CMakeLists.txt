add_subdirectory(minimizer)
add_subdirectory(models)
if(Boost_FOUND)
  add_subdirectory(baryo_calculation)
  add_subdirectory(WallThickness)
  add_subdirectory(Kfactors)
endif(Boost_FOUND)
add_subdirectory(ThermalFunctions)


#==============================Utility==================================

set(header_path "${CMAKE_SOURCE_DIR}/include/BSMPT")
set(header ${CMAKE_SOURCE_DIR}/include/BSMPT/utility.h)
set(src utility.cpp)
add_library(Utility ${header} ${src})
target_include_directories(Utility PUBLIC ${BSMPT_SOURCE_DIR}/include)


#==============================EXE======================================

add_executable(BSMPT prog/BSMPT.cpp)
target_link_libraries(BSMPT Minimizer Models Utility)
target_compile_features(BSMPT PUBLIC cxx_std_14)

add_executable(VEVEVO prog/VEVEVO.cpp)
target_link_libraries(VEVEVO Minimizer Models Utility)
target_compile_features(VEVEVO PUBLIC cxx_std_14)

add_executable(NLOVEV prog/NLOVEV.cpp)
target_link_libraries(NLOVEV Minimizer Models Utility)
target_compile_features(NLOVEV PUBLIC cxx_std_14)

add_executable(CalcCT prog/CalcCT.cpp)
target_link_libraries(CalcCT Models Utility)
target_compile_features(CalcCT PUBLIC cxx_std_14)

add_executable(TripleHiggsCouplingsNLO prog/TripleHiggsNLO.cpp)
target_link_libraries(TripleHiggsCouplingsNLO Models Utility)
target_compile_features(TripleHiggsCouplingsNLO PUBLIC cxx_std_14)

add_executable(Test prog/Test.cpp)
target_link_libraries(Test Minimizer Models Utility)
target_compile_features(Test PUBLIC cxx_std_14)



if(Boost_FOUND)
add_executable(CalculateEWBG prog/CalculateEWBG.cpp)
target_link_libraries(CalculateEWBG Minimizer Models Baryo Utility)
target_compile_features(CalculateEWBG PUBLIC cxx_std_14)
else(Boost_FOUND)
  message("CalculateEWBG will not compile as Boost is missing")
endif(Boost_FOUND)

if(Boost_FOUND)
add_executable(PlotEWBG_vw prog/PlotEWBG_vw.cpp)
target_link_libraries(PlotEWBG_vw Minimizer Models Baryo Utility)
target_compile_features(PlotEWBG_vw PUBLIC cxx_std_14)
else(Boost_FOUND)
  message("PlotEWBG_vw will not compile as Boost is missing")
endif(Boost_FOUND)

if(Boost_FOUND)
add_executable(PlotEWBG_nL prog/PlotEWBG_nL.cpp)
target_link_libraries(PlotEWBG_nL Minimizer Models Baryo Utility)
target_compile_features(PlotEWBG_nL PUBLIC cxx_std_14)
else(Boost_FOUND)
  message("PlotEWBG_nL will not compile as Boost is missing")
endif(Boost_FOUND)

add_executable(RenormScale prog/RenormScale.cpp)
target_link_libraries(RenormScale Minimizer Models Utility)
target_compile_features(RenormScale PUBLIC cxx_std_14)

if(Boost_FOUND)
add_executable(EWBGRenormScale prog/PlotEWBG_mu.cpp)
target_link_libraries(EWBGRenormScale Minimizer Models Baryo Utility)
target_compile_features(EWBGRenormScale PUBLIC cxx_std_14)
else(Boost_FOUND)
  message("EWBGRenormScale will not compile as Boost is missing")
endif(Boost_FOUND)



#======================= Additional testing ====================================

if(MakeAdditionalTesting)

  if(Boost_FOUND)
    add_executable(CreateMuGrid prog/CreateMuGrid.cpp)
    target_link_libraries(CreateMuGrid Minimizer Models Baryo Utility)
    target_compile_features(CreateMuGrid PUBLIC cxx_std_14)
  else(Boost_FOUND)
    message("CreateMuGrid will not compile as Boost is missing")
  endif(Boost_FOUND)

  if(Boost_FOUND)
    add_executable(Kfactor_grid prog/Kfactor_grid.cpp)
    target_link_libraries(Kfactor_grid Minimizer Models Kfactors Utility)
    target_compile_features(Kfactor_grid PUBLIC cxx_std_14)
  else(Boost_FOUND)
    message("Kfactor_grid will not compile as Boost is missing")
  endif(Boost_FOUND)

  if(Boost_FOUND)
  add_executable(WallThickness prog/WallThickness.cpp)
  target_link_libraries(WallThickness Minimizer Models LibWallThickness Utility)
  target_compile_features(WallThickness PUBLIC cxx_std_14)
  else(Boost_FOUND)
    message("WallThickness will not compile as Boost is missing")
  endif(Boost_FOUND)

endif(MakeAdditionalTesting)

