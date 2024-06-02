// Copyright (C) 2018  Philipp Basler and Margarete Mühlleitner
// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include <BSMPT/models/ClassPotentialC2HDM.h>
#include <BSMPT/models/ClassPotentialCPintheDark.h>
#include <BSMPT/models/ClassPotentialCxSM.h>
#include <BSMPT/models/ClassPotentialN2HDM.h>
#include <BSMPT/models/ClassPotentialOrigin.h> // for Class_Potential_Origin
#include <BSMPT/models/ClassPotentialR2HDM.h>
#include <BSMPT/models/ClassPotentialSM.h>
#include <BSMPT/models/IncludeAllModels.h>
#include <ctype.h>   // for isdigit, tolower
#include <iostream>  // for operator<<, cerr, ost...
#include <stdexcept> // for runtime_error
#include <utility>   // for pair

#include <BSMPT/models/ClassTemplate.h>
#include <BSMPT/utility/Logger.h>

namespace BSMPT
{
namespace ModelID
{

std::unique_ptr<Class_Potential_Origin> FChoose(ModelIDs choice)
{
  return FChoose(choice, GetSMConstants());
}

std::unique_ptr<Class_Potential_Origin> FChoose(ModelIDs choice,
                                                const ISMConstants &smConstants)
{
  using namespace Models;
  switch (choice)
  {
  case ModelIDs::SM: return std::make_unique<Class_SM>(smConstants); break;
  case ModelIDs::R2HDM:
    return std::make_unique<Class_Potential_R2HDM>(smConstants);
    break;
  case ModelIDs::C2HDM:
    return std::make_unique<Class_Potential_C2HDM>(smConstants);
    break;
  case ModelIDs::N2HDM:
    return std::make_unique<Class_Potential_N2HDM>(smConstants);
    break;
  case ModelIDs::CXSM: return std::make_unique<Class_CxSM>(smConstants); break;
  case ModelIDs::CPINTHEDARK:
    return std::make_unique<Class_Potential_CPintheDark>(smConstants);
    break;
  case ModelIDs::TEMPLATE:
    return std::make_unique<Class_Template>(smConstants);
    break;
  default: throw std::runtime_error("Invalid model");
  }
}

ModelIDs getModel(const std::string &s)
{
  std::string ModelInput = s;
  std::transform(
      ModelInput.begin(), ModelInput.end(), ModelInput.begin(), ::tolower);

  for (const auto &entry : ModelNames)
  {
    auto key = entry.first;
    std::transform(key.begin(), key.end(), key.begin(), ::tolower);
    if (ModelInput.compare(key) == 0)
    {
      return entry.second;
    }
  }
  return ModelIDs::NotSet;
}

} // namespace ModelID

} // namespace BSMPT
