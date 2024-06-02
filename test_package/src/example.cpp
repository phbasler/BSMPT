#include <BSMPT/utility/utility.h>
#include <BSMPT/utility/asciiplotter/asciiplotter.h>
#include <BSMPT/utility/spline/spline.h>
#include <BSMPT/bounce_solution/action_calculation.h>
#include <BSMPT/gravitational_waves/gw.h>
#include <BSMPT/minimizer/Minimizer.h>
#include <BSMPT/minimum_tracer/minimum_tracer.h>
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/ThermalFunctions/ThermalFunctions.h>
#include <BSMPT/transition_tracer/transition_tracer.h>

int main() {
    using namespace BSMPT;
    Logger::Disable();
    ModelID::FChoose(ModelID::ModelIDs::C2HDM, GetSMConstants());
    Minimizer::CalcWhichMinimizer();

}
