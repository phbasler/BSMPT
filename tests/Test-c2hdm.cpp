#include "catch.hpp"

#include <BSMPT/models/ClassPotentialOrigin.h>  // for Class_Potential_Origin
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/minimizer/Minimizer.h>

const std::vector<double> example_point_C2HDM{
    /* lambda_1 = */ 3.29771,
    /* lambda_2 = */ 0.274365,
    /* lambda_3 = */ 4.71019,
    /* lambda_4 = */ -2.23056,
    /* Re(lambda_5) = */ -2.43487,
    /* Im(lambda_5) = */ 0.124948,
    /* Re(m_{12}^2) = */ 2706.86,
    /* tan(beta) = */ 4.64487,
    /* Yukawa Type = */ 1
};

TEST_CASE("Checking NLOVEV for C2HDM", "[c2hdm]") {
    using namespace BSMPT;
    std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer = ModelID::FChoose(ModelID::ModelIDs::C2HDM);
    modelPointer->initModel(example_point_C2HDM);
    std::vector<double> Check;
    auto sol = Minimizer::Minimize_gen_all(modelPointer,0,Check,modelPointer->get_vevTreeMin(),Minimizer::WhichMinimizerDefault);
    for(std::size_t i{0};i<sol.size();++i)
    {
        auto expected = std::abs(modelPointer->get_vevTreeMin(i));
        auto res = std::abs(sol.at(i));
        REQUIRE( std::abs(res-expected) <= 1e-4);
    }
}


TEST_CASE("Checking EWPT for C2HDM", "[c2hdm]") {
    using namespace BSMPT;
    std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer = ModelID::FChoose(ModelID::ModelIDs::C2HDM);
    modelPointer->initModel(example_point_C2HDM);
    std::vector<double> Check;
    auto EWPT = Minimizer::PTFinder_gen_all(modelPointer,0,300,Minimizer::WhichMinimizerDefault);
    const double omega_c_expected = 200.79640966130026;
    const double Tc_expected = 145.56884765625;
    const std::vector<double> min_expected{0, -49.929666284908336, -194.48507400070247, 1.3351211509361505};
    REQUIRE(EWPT.StatusFlag == Minimizer::MinimizerStatus::SUCCESS);
    REQUIRE(std::abs(omega_c_expected - EWPT.vc)/omega_c_expected <= 1e-4);
    REQUIRE(std::abs(Tc_expected-EWPT.Tc)/Tc_expected <= 1e-4);
    for(std::size_t i{0};i<EWPT.EWMinimum.size();++i)
    {
        auto res = std::abs(EWPT.EWMinimum.at(i));
        auto expected = std::abs(min_expected.at(i));
        if(expected != 0)
        {
            REQUIRE(std::abs(res-expected)/expected <= 1e-4);
        }
        else {
            REQUIRE(res <= 1e-4);
        }
    }
}
