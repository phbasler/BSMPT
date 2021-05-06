// Copyright (C) 2020  Philipp Basler, Margarete Mühlleitner and Jonas Müller
// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/utility.h>
#include <string>

/**
 * @file
 */
namespace BSMPT
{

std::ostream &operator<<(std::ostream &os, const ModelID::ModelIDs &Model)
{
  static auto IMN = BSMPT::ModelID::InvertModelNames();
  os << IMN.at(Model);
  return os;
}
bool StringStartsWith(const std::string &str, const std::string &prefix)
{
  return str.size() >= prefix.size() and str.find(prefix) == 0;
}

} // namespace BSMPT
