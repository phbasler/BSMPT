// Copyright (C) 2020  Philipp Basler, Margarete Mühlleitner and Jonas Müller
// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

/**
 * @file
 */

#include <iostream>

#include <BSMPT/utility/Logger.h>
namespace BSMPT
{
void BSMPTLogger::SetOStream(std::ostream &Ostream)
{
  mfilestream.close();
  mOstream.rdbuf(Ostream.rdbuf());
}

void BSMPTLogger::SetLevel(const std::map<LoggingLevel, bool> &Setup)
{
  for (const auto &el : Setup)
  {
    auto pos = mCurrentSetup.find(el.first);
    if (pos != mCurrentSetup.end())
    {
      pos->second = el.second;
    }
    else
    {
      mCurrentSetup.emplace(el);
    }
  }
}

void BSMPTLogger::SetOStream(const std::string &file)
{
  mfilestream = std::ofstream(file);
  mOstream.rdbuf(mfilestream.rdbuf());
}

void BSMPTLogger::Disable()
{
  mCurrentSetup.clear();
}

} // namespace BSMPT
