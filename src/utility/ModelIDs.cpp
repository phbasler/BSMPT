#include <BSMPT/utility/Logger.h>
#include <BSMPT/utility/ModelIDs.h>
#include <BSMPT/utility/parser.h>
#include <BSMPT/utility/settings.h>
#include <sstream>
#include <stdexcept> // for runtime_error

namespace BSMPT
{
namespace ModelID
{

std::unordered_map<ModelIDs, std::string> InvertModelNames()
{
  std::unordered_map<ModelIDs, std::string> IMN;
  for (const auto &el : ModelNames)
  {
    auto success = IMN.emplace(el.second, el.first);
    if (not success.second)
    {
      throw std::runtime_error(
          "\nERROR: The same ModelID is assigned for two different models.\n");
    }
  }
  return IMN;
}

} // namespace ModelID



std::ostream &operator<<(std::ostream &os, const ModelID::ModelIDs &Model)
{
  static auto IMN = BSMPT::ModelID::InvertModelNames();
  os << IMN.at(Model);
  return os;
}

} // namespace BSMPT