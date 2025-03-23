// Copyright (C) 2020  Philipp Basler, Margarete Mühlleitner and Jonas Müller
// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include <fstream>
#include <iostream>
#include <map>
#include <vector>

/**
 * @file
 */
namespace BSMPT
{

class Logger;
class parser;

void ShowLoggerHelp();
void SetLogger(const parser &argparser);
void SetLogger(const std::vector<std::string> &args);

enum class LoggingLevel
{
  None,
  Default,
  MinimizerDetailed,
  MinTracerDetailed,
  TransitionDetailed,
  GWDetailed,
  BounceDetailed,
  ProgDetailed,
  EWBGDetailed,
  Debug,
  Complete
};

/**
 * @brief The BSMPTLogger class
 */
class BSMPTLogger
{
public:
  friend Logger;
  BSMPTLogger(std::ostream &os) : mOstream{os.rdbuf()} {}
  BSMPTLogger(const BSMPTLogger &)            = delete;
  BSMPTLogger(BSMPTLogger &&)                 = delete;
  BSMPTLogger &operator=(const BSMPTLogger &) = delete;
  BSMPTLogger &operator=(BSMPTLogger &&)      = delete;

private:
  /**
   * @brief SetOStream sets the output of the logger to a certain ostream
   * @param Ostream
   */
  void SetOStream(std::ostream &Ostream);
  /**
   * @brief SetOStream writes the output of the logger to a file
   * @param file
   */
  void SetOStream(const std::string &file);
  template <typename T>
  void
  Write(LoggingLevel level, const T &toWrite, const std::string &file, int line)
  {
    auto pos = mCurrentSetup.find(level);
    if (pos != mCurrentSetup.end() and pos->second)
    {
      if (not file.empty())
      {
        mOstream << "file: " << file << "; ";
      }
      if (line >= 0)
      {
        mOstream << "Line: " << line << "; ";
      }
      mOstream << toWrite << std::endl;
    }
  }
  void SetLevel(const std::map<LoggingLevel, bool> &level);
  void SetLevel(LoggingLevel level, bool enable);
  void Disable();
  void RestoreDefaultLevels();

  std::ostream mOstream;
  std::ofstream mfilestream;

  const std::map<LoggingLevel, bool> mDefaultSetup{
      {LoggingLevel::Default, true},
      {LoggingLevel::EWBGDetailed, false},
      {LoggingLevel::ProgDetailed, false},
      {LoggingLevel::MinimizerDetailed, false},
      {LoggingLevel::MinTracerDetailed, false},
      {LoggingLevel::TransitionDetailed, false},
      {LoggingLevel::BounceDetailed, false},
      {LoggingLevel::GWDetailed, false},
      {LoggingLevel::Debug, false}};

  std::map<LoggingLevel, bool> mCurrentSetup = mDefaultSetup;
};

class Logger
{
public:
  Logger(const Logger &)            = delete;
  Logger(Logger &&)                 = delete;
  Logger &operator=(const Logger &) = delete;
  Logger &operator=(Logger &&)      = delete;

  static void SetLevel(const std::map<LoggingLevel, bool> &Setup)
  {
    Instance().SetLevel(Setup);
  }
  static void SetLevel(LoggingLevel level, bool enable)
  {
    Instance().SetLevel(level, enable);
  }

  static void SetOStream(std::ostream &Ostream)
  {
    Instance().SetOStream(Ostream);
  }
  static void SetOStream(const std::string &file)
  {
    Instance().SetOStream(file);
  }
  template <typename T>
  static void Write(LoggingLevel level,
                    const T &toWrite,
                    const std::string &file = "",
                    int line                = -1)
  {
    Instance().Write(level, toWrite, file, line);
  }

  static bool GetLoggingLevelStatus(LoggingLevel level)
  {
    return Instance().mCurrentSetup[level];
  }

  static void Disable() { Instance().Disable(); }

  static void RestoreDefaultLevels() { Instance().RestoreDefaultLevels(); }

private:
  static BSMPTLogger &Instance()
  {
    static BSMPTLogger LoggerInstance(std::cout);
    return LoggerInstance;
  }
};

} // namespace BSMPT
