#include <BSMPT/ThermalFunctions/ThermalFunctions.h>
#include <BSMPT/bounce_solution/action_calculation.h>
#include <BSMPT/gravitational_waves/gw.h>
#include <BSMPT/minimizer/Minimizer.h>
#include <BSMPT/minimum_tracer/minimum_tracer.h>
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/transition_tracer/transition_tracer.h>
#include <BSMPT/utility/asciiplotter/asciiplotter.h>
#include <BSMPT/utility/spline/spline.h>
#include <BSMPT/utility/utility.h>

int main()
{
  using namespace BSMPT;
  Logger::Disable();
  std::shared_ptr<Class_Potential_Origin> testModel =
      ModelID::FChoose(ModelID::ModelIDs::C2HDM, GetSMConstants());
  Minimizer::CalcWhichMinimizer();

  const std::vector<double> example_point_C2HDM{/* lambda_1 = */ 3.29771,
                                                /* lambda_2 = */ 0.274365,
                                                /* lambda_3 = */ 4.71019,
                                                /* lambda_4 = */ -2.23056,
                                                /* Re(lambda_5) = */ -2.43487,
                                                /* Im(lambda_5) = */ 0.124948,
                                                /* Re(m_{12}^2) = */ 2706.86,
                                                /* tan(beta) = */ 4.64487,
                                                /* Yukawa Type = */ 1};

  testModel->initModel(example_point_C2HDM);

  std::vector<double> check;
  std::vector<double> start = testModel->get_vevTreeMin();
  auto result =
      Minimizer::Minimize_gen_all(testModel,
                                  0,
                                  check,
                                  start,
                                  BSMPT::Minimizer::WhichMinimizerDefault,
                                  false);
  std::cout << "result = " << result << std::endl;
}
