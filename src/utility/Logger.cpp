// Copyright (C) 2020  Philipp Basler, Margarete Mühlleitner and Jonas Müller
// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

/**
 * @file
 */

#include <BSMPT/utility/Logger.h>
#include <BSMPT/utility/parser.h>
#include <BSMPT/utility/utility.h>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <sstream>

namespace BSMPT
{

static std::map<std::string, LoggingLevel> LoggingPrefixes{
    {"--logginglevel::default=", LoggingLevel::Default},
    {"--logginglevel::debug=", LoggingLevel::Debug},
    {"--logginglevel::ewbgdetailed=", LoggingLevel::EWBGDetailed},
    {"--logginglevel::progdetailed=", LoggingLevel::ProgDetailed},
    {"--logginglevel::minimizerdetailed=", LoggingLevel::MinimizerDetailed},
    {"--logginglevel::transitiondetailed=", LoggingLevel::TransitionDetailed},
    {"--logginglevel::mintracerdetailed=", LoggingLevel::MinTracerDetailed},
    {"--logginglevel::bouncedetailed=", LoggingLevel::BounceDetailed},
    {"--logginglevel::gwdetailed=", LoggingLevel::GWDetailed},
    {"--logginglevel::complete=", LoggingLevel::Complete}};

void ShowLoggerHelp()
{
  int SizeOfFirstColumn =
      std::string("--logginglevel::minimizerdetailed=           ").size();
  std::stringstream ss;
  ss << std::setw(SizeOfFirstColumn) << std::left
     << "The following options for the Logger are available:" << std::endl
     << std::setw(SizeOfFirstColumn) << std::left
     << "--logginglevel::disabled to disable the Logger." << std::endl;
  for (const auto &el : LoggingPrefixes)
  {
    ss << std::setw(SizeOfFirstColumn) << std::left << el.first
       << "\t true/false" << std::endl;
  }
  Logger::Write(LoggingLevel::Default, ss.str());
}

void SetLogger(const BSMPT::parser &argparser)
{
  try
  {
    argparser.get_value("logginglevel::disabled");
    Logger::Disable();
    return;
  }
  catch (BSMPT::parserException &)
  {
  }

  for (const auto &[prefix, enumValue] : LoggingPrefixes)
  {
    auto beginning  = prefix.find_first_not_of("-");
    auto seperator  = prefix.find("=");
    std::string key = prefix.substr(beginning, seperator - beginning);

    try
    {
      auto value = argparser.get_value(key);
      Logger::SetLevel(enumValue, value == "true");
      // Catch logginglevel::complete
      if ((key == "logginglevel::complete") and (value == "true"))
      {
        for (const auto &[prefix_, enumValue_] : LoggingPrefixes)
          Logger::SetLevel(enumValue_, true);
        Logger::SetLevel(LoggingLevel::MinimizerDetailed, false);
      }
    }
    catch (BSMPT::parserException &)
    {
    }
  }
}

void SetLogger(const std::vector<std::string> &args)
{
  auto posDisable =
      std::find(args.begin(), args.end(), "--logginglevel::disabled");
  if (posDisable != args.end())
  {
    Logger::Disable();
    return;
  }

  for (const auto &el : args)
  {
    auto pos =
        std::find_if(LoggingPrefixes.begin(),
                     LoggingPrefixes.end(),
                     [&el](auto pr) { return StringStartsWith(el, pr.first); });
    if (pos != LoggingPrefixes.end())
    {
      Logger::SetLevel(pos->second, el.substr(pos->first.size()) == "true");
    }
    // Catch logginglevel::complete
    if ((pos->second == LoggingLevel::Complete) and
        (el.substr(pos->first.size()) == "true"))
    {
      for (const auto &[prefix, enumValue] : LoggingPrefixes)
        Logger::SetLevel(enumValue, true);
      Logger::SetLevel(LoggingLevel::MinimizerDetailed, false);
    }
  }
}

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

void BSMPTLogger::SetLevel(LoggingLevel level, bool enable)
{
  mCurrentSetup[level] = enable;
}

void BSMPTLogger::SetOStream(const std::string &file)
{
  mfilestream = std::ofstream(file);
  mOstream.rdbuf(mfilestream.rdbuf());
}

void BSMPTLogger::Disable()
{
  std::cout << "Disable logger" << std::endl;
  mCurrentSetup.clear();
}

void BSMPTLogger::RestoreDefaultLevels()
{
  mCurrentSetup = mDefaultSetup;
}

} // namespace BSMPT
