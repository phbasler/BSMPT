#pragma once

#include <string>
#include <unordered_map>

namespace BSMPT
{
namespace ModelID
{

/**
 * @brief The ModelIDs enum containing all IDs for identifying the Models
 */
enum class ModelIDs
{
  NotSet,
  C2HDM,
  R2HDM,
  N2HDM,
  CXSM,
  CPINTHEDARK,
  SM,

  // Here you start adding your models
  TEMPLATE,

  // DO NOT EDIT the part below
  stop
};

/**
 * @brief Mapping between the model name which is given as the first argument to
 * the binary and the ModelIDs element
 */
const std::unordered_map<std::string, ModelIDs> ModelNames{
    {"c2hdm", ModelIDs::C2HDM},
    {"r2hdm", ModelIDs::R2HDM},
    {"n2hdm", ModelIDs::N2HDM},
    {"cxsm", ModelIDs::CXSM},
    {"sm", ModelIDs::SM},
    {"cpinthedark", ModelIDs::CPINTHEDARK},
    {"template", ModelIDs::TEMPLATE},
};

} // namespace ModelID

/**
 * @brief operator << overload for the model parameter
 */
std::ostream &operator<<(std::ostream &os, const ModelID::ModelIDs &Model);

} // namespace BSMPT