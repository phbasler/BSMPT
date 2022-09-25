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
  ProgDetailed,
  EWBGDetailed,
  Debug
};

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
  void SetOStream(std::ostream &Ostream);
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

  std::ostream mOstream;
  std::ofstream mfilestream;

  std::map<LoggingLevel, bool> mCurrentSetup{
      {LoggingLevel::Default, true},
      {LoggingLevel::EWBGDetailed, false},
      {LoggingLevel::ProgDetailed, false},
      {LoggingLevel::MinimizerDetailed, false},
      {LoggingLevel::Debug, false}};
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

  static void Disable() { Instance().Disable(); }

private:
  static BSMPTLogger &Instance()
  {
    static BSMPTLogger LoggerInstance(std::cout);
    return LoggerInstance;
  }
};

} // namespace BSMPT
