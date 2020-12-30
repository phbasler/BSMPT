#include "catch.hpp"

#include <BSMPT/models/ClassPotentialOrigin.h>  // for Class_Potential_Origin
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/minimizer/Minimizer.h>
#include <BSMPT/utility.h>

const std::vector<double> example_point_RN2HDM{
    /* lambda_1 = */ 0.300812 ,
    /* lambda_2 = */ 0.321809,
    /* lambda_3 = */ -0.133425,
    /* lambda_4 = */ 4.11105,
    /* lambda_5 = */ -3.84178,
    /* lambda_6 = */ 9.46329,
    /* lambda_7 = */ -0.750455,
    /* lambda_8 = */ 0.743982,
    /* tan(beta) = */ 5.91129,
    /* v_s = */ 293.035,
    /* m_{12}^2 = */ 4842.28,
    /* Yukawa Type = */ 1
};

constexpr auto Model = BSMPT::ModelID::ModelIDs::RN2HDM;

TEST_CASE("Checking NLOVEV for N2HDM", "[n2hdm]") {
    using namespace BSMPT;
    std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer = ModelID::FChoose(Model);
    modelPointer->initModel(example_point_RN2HDM);
    std::vector<double> Check;
    auto sol = Minimizer::Minimize_gen_all(modelPointer,0,Check,modelPointer->get_vevTreeMin(),Minimizer::WhichMinimizerDefault);
    for(std::size_t i{0};i<sol.size();++i)
    {
        auto expected = std::abs(modelPointer->get_vevTreeMin(i));
        auto res = std::abs(sol.at(i));
        REQUIRE( std::abs(res-expected) <= 1e-4);
    }
}


TEST_CASE("Checking EWPT for N2HDM", "[N2hdm]") {
    using namespace BSMPT;
    std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer = ModelID::FChoose(Model);
    modelPointer->initModel(example_point_RN2HDM);
    std::vector<double> Check;
    auto EWPT = Minimizer::PTFinder_gen_all(modelPointer,0,300,Minimizer::WhichMinimizerDefault);
    const double omega_c_expected = 180.5917335676395;
    const double Tc_expected = 120.7305908203125;
    const std::vector<double> min_expected{0, 0, -32.70827526931041, -177.6050195289305, -297.0418903961274};
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
