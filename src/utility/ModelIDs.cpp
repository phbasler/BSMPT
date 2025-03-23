#include <BSMPT/utility/Logger.h>
#include <BSMPT/utility/ModelIDs.h>
#include <BSMPT/utility/parser.h>
#include <BSMPT/utility/settings.h>
#include <BSMPT/utility/utility.h>
#include <sstream>
#include <stdexcept> // for runtime_error

namespace BSMPT
{

std::ostream &operator<<(std::ostream &os, const ModelID::ModelIDs &Model)
{
  static auto IMN = BSMPT::InvertMap(
      ModelID::ModelNames,
      "\nERROR: The same ModelID is assigned for two different models.\n");
  os << IMN.at(Model);
  return os;
}

} // namespace BSMPT