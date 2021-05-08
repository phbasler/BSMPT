// Copyright (C) 2020  Philipp Basler, Margarete Mühlleitner and Jonas Müller
// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include <fstream>

/**
 * @file
 */
namespace BSMPT
{
class Logger;

enum class LoggingLevel
{
  None,
  Default,
  Detailed
};

class BSMPTLogger
{
public:
  friend Logger;
  BSMPTLogger(std::ostream &os) : mOstream{os.rdbuf()} {}

private:
  BSMPTLogger(const BSMPTLogger &) = default;
  BSMPTLogger(BSMPTLogger &&)      = default;
  BSMPTLogger &operator=(const BSMPTLogger &) = default;
  BSMPTLogger &operator=(BSMPTLogger &&) = default;

  void SetOStream(std::ostream &Ostream);
  template <typename T> void Write(LoggingLevel level, const T &toWrite)
  {
    if (level <= mCurrentLevel)
    {
      mOstream << toWrite << std::endl;
    }
  }
  void SetLevel(LoggingLevel level);
  void SetLoggingFile(const std::string &file);

  std::ostream mOstream;
  std::ofstream mfilestream;
  LoggingLevel mCurrentLevel{LoggingLevel::Default};
};

class Logger
{
public:
  Logger(const Logger &) = delete;
  Logger(Logger &&)      = delete;
  Logger &operator=(const Logger &) = delete;
  Logger &operator=(Logger &&) = delete;

  static void SetLevel(LoggingLevel level) { Instance().SetLevel(level); }

  static void SetOStream(std::ostream &Ostream)
  {
    Instance().SetOStream(Ostream);
  }
  template <typename T> static void Write(LoggingLevel level, const T &toWrite)
  {
    Instance().Write(level, toWrite);
  }

  static void SetLoggingFile(const std::string &file)
  {
    Instance().SetLoggingFile(file);
  }

private:
  static BSMPTLogger &Instance()
  {
    static BSMPTLogger LoggerInstance(std::cout);
    return LoggerInstance;
  }
};

} // namespace BSMPT
