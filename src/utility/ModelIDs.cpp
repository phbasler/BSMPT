#include <BSMPT/utility/ModelIDs.h>
#include <stdexcept> // for runtime_error
#include <sstream>
#include <BSMPT/utility/Logger.h>
#include <BSMPT/utility/parser.h>
#include <BSMPT/utility/settings.h>

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

void ShowInputError()
{
  std::stringstream ss;
  ss << "The chosen Method for the thermal mass corrections is ";
  if (C_UseParwani)
    ss << "Parwani ";
  else
    ss << "Arnold Espinosa\n";
  ss << "The implemented models are " << std::endl;
  for (auto entry : ModelID::ModelNames)
  {
    ss << entry.first << std::endl;
  }
  Logger::Write(LoggingLevel::Default, ss.str());
}

std::ostream &operator<<(std::ostream &os, const ModelID::ModelIDs &Model)
{
  static auto IMN = BSMPT::ModelID::InvertModelNames();
  os << IMN.at(Model);
  return os;
}

std::string ModelIDToString(const ModelID::ModelIDs &Model)
{
  std::stringstream ss;
  ss << Model;
  return ss.str();
}

} // namespace BSMPT