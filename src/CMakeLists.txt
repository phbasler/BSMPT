add_subdirectory(minimizer)
add_subdirectory(models)

add_executable(BSMPT prog/BSMPT.cpp)
target_link_libraries(BSMPT Minimizer Models)

add_executable(VEVEVO prog/VEVEVO.cpp)
target_link_libraries(VEVEVO Minimizer Models)

add_executable(NLOVEV prog/NLOVEV.cpp)
target_link_libraries(NLOVEV Minimizer Models)

add_executable(CalcCT prog/CalcCT.cpp)
target_link_libraries(CalcCT Models)

add_executable(TripleHiggsCouplingsNLO prog/TripleHiggsNLO.cpp)
target_link_libraries(TripleHiggsCouplingsNLO Models)

add_executable(Test prog/Test.cpp)
target_link_libraries(Test Minimizer Models)
