// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include <BSMPT/config.h>
#include <BSMPT/utility/Logger.h>
#include <BSMPT/utility/ModelIDs.h>
#include <BSMPT/utility/parser.h>
#include <BSMPT/utility/settings.h>
#include <BSMPT/utility/utility.h>
#include <algorithm>
#include <sstream>

#ifdef nlohmann_json_FOUND
#include <fstream>
#include <nlohmann/json.hpp>
#endif

namespace BSMPT
{

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

parserException::parserException(const std::string &msg) : message{msg}
{
}

const char *parserException::what() const noexcept
{
  return message.c_str();
}

parser::parser()
{
  add_argument("logginglevel::default",
               "Turn on/off (true/false) the default output.",
               false);
  add_argument("logginglevel::debug",
               "Turn on/off (true/false) the debug output for the logger.",
               false);
  add_argument("logginglevel::disabled", "Disable the Logger.", false);

  add_argument("logginglevel::ewbgdetailed",
               "Turn on/off (true/false) the output for the EWBG calculation.",
               false);
  add_argument(
      "logginglevel::progdetailed",
      "Turn on/off (true/false) the additional output for the main program.",
      false);
  add_argument(
      "logginglevel::minimizerdetailed",
      "Turn on/off (true/false) the additional output of the minimizers.",
      false);

  add_argument("json",
               "Use a json file to define the input instead of additional cli "
               "parameters.",
               false);
}

parser::parser(const bool &enable_column_output)
{
  this->extra_column_output = enable_column_output;

  add_argument("logginglevel::default", "default output", "true", false);
  add_argument("logginglevel::debug", "debug output", "false", false);
  add_argument("logginglevel::disabled", "disable all output", "", false);
  add_argument("logginglevel::ewbgdetailed",
               "baryogenesis calculation output",
               "false",
               false);
  add_argument(
      "logginglevel::progdetailed", "executable output", "false", false);
  add_argument(
      "logginglevel::minimizerdetailed", "minimizer output", "false", false);
  add_argument("logginglevel::transitiondetailed",
               "transition calculation output",
               "false",
               false);
  add_argument("logginglevel::mintracerdetailed",
               "minimum tracking output",
               "false",
               false);
  add_argument("logginglevel::bouncedetailed",
               "bounce solution calculation output",
               "false",
               false);
  add_argument("logginglevel::gwdetailed",
               "gravitational wave calculation output",
               "false",
               false);
  add_argument("logginglevel::complete",
               "enable all except minimizer output",
               "false",
               false);
}

void parser::enable_minimizer_options()
{
  add_argument("useGSL",
               "Use the GSL library to minimize the effective potential.",
               false);
  add_argument("useCMAES",
               "Use the CMAES library to minimize the effective potential",
               false);
  add_argument("useNLopt",
               "Use the NLopt library to minimize the effective potential",
               false);
  add_argument("useMultithreading",
               "Enables/Disables multi threading for the minimizers",
               false);
}

void parser::add_argument_only_display(const std::string &argument,
                                       const std::string &description,
                                       const std::string &default_val)
{
  auto arg = to_lower(argument);
  Options options;
  options.argument    = argument;
  options.description = description;
  options.default_val = default_val;

  mOrderedArguments.push_back(std::pair<std::string, Options>(arg, options));
}

void parser::add_argument(const std::string &argument, bool required)
{
  auto arg = to_lower(argument);
  Options options;
  options.argument = argument;
  options.value    = std::nullopt;

  if (required)
  {
    mRequiredArguments.emplace(arg, options);
  }
  else
  {
    mOptionalArguments.emplace(arg, options);
  }
}

void parser::add_argument(const std::string &argument,
                          const std::string &description,
                          bool required)
{
  auto arg = to_lower(argument);
  Options options;
  options.argument    = argument;
  options.description = description;
  options.value       = std::nullopt;

  mOrderedArguments.push_back(std::pair<std::string, Options>(arg, options));

  if (required)
  {
    mRequiredArguments.emplace(arg, options);
  }
  else
  {
    mOptionalArguments.emplace(arg, options);
  }
}

void parser::add_argument(const std::string &argument,
                          const std::string &description,
                          const std::string &default_val,
                          bool required)
{
  auto arg = to_lower(argument);
  Options options;
  options.argument    = argument;
  options.description = description;
  options.value       = std::nullopt;
  options.default_val = default_val;

  mOrderedArguments.push_back(std::pair<std::string, Options>(arg, options));

  if (required)
  {
    mRequiredArguments.emplace(arg, options);
  }
  else
  {
    mOptionalArguments.emplace(arg, options);
  }
}

void parser::add_subtext(const std::string &subtext)
{
  Options options;
  options.description = subtext;
  mOrderedArguments.push_back(
      std::pair<std::string, Options>("subtext", options));
}

void parser::add_input(const std::vector<std::string> &input)
{
  if (input.size() == 1)
  {
    std::string arg{input.at(0)};
    auto [argument, value] = get_key_value(arg);
    if (argument == "help")
    {
      print_help();
      return;
    }
    else if (argument == "json")
    {
#ifdef nlohmann_json_FOUND
      add_json_input(value);
      return;
#else
      throw parserException(
          "nlohmann_json is required for the json config file.");
#endif
    }
  }

  std::vector<KeyValue> sepInput;
  for (const auto &arg : input)
  {
    sepInput.emplace_back(get_key_value(arg));
  }
  add_input(sepInput);
}

bool parser::all_required_set() const
{
  for (const auto &[key, value] : mRequiredArguments)
  {
    (void)key;
    if (not value.value.has_value())
    {
      Logger::Write(LoggingLevel::Default,
                    "The required parameter " + key + " is not set.");
      return false;
    }
  }
  return true;
}

void parser::check_required_parameters() const
{
  if (not all_required_set())
  {
    print_help();
    throw parserException("Not all required parameters are set.");
  }
}

void parser::print_help() const
{
  if (mHelpAlreadyPrinted)
  {
    return;
  }
  mHelpAlreadyPrinted = true;
  std::stringstream ss, ss_end;
  ss << mHeader << std::endl;
  ss_end << "The following options for the Logger are available:\n\n";

  int size_first_column  = 37;
  int size_second_column = 10;
  if (extra_column_output)
  {
    ss << std::setw(size_first_column) << std::left << "argument"
       << std::setw(size_second_column) << std::left << "default"
       << "description" << std::endl;

    for (const auto &el : mOrderedArguments)
    {
      if (el.first == "subtext")
      {
        ss << std::setw(size_first_column) << std::left << ""
           << std::setw(size_second_column) << std::left << ""
           << el.second.description << std::endl;
      }
      else if (StringStartsWith(el.second.argument, "logginglevel"))
      {
        if (el.second.argument.find("disable"))
        {
          ss_end << std::setw(size_first_column) << std::left
                 << "--" + el.second.argument << std::setw(size_second_column)
                 << std::left << el.second.default_val << el.second.description
                 << std::endl;
        }
        else
        {
          ss_end << std::setw(size_first_column) << std::left
                 << "--" + el.second.argument + "="
                 << std::setw(size_second_column) << std::left
                 << el.second.default_val << el.second.description << std::endl;
        }
      }
      else if (StringStartsWith(el.second.argument, "help"))
      {
        ss << std::setw(size_first_column) << std::left
           << "--" + el.second.argument << std::setw(size_second_column)
           << std::left << el.second.default_val << el.second.description
           << std::endl;
      }
      else
      {
        ss << std::setw(size_first_column) << std::left
           << "--" + el.second.argument + "=" << std::setw(size_second_column)
           << std::left << el.second.default_val << el.second.description
           << std::endl;
      }
    }
    Logger::Write(LoggingLevel::Default, ss.str());
    Logger::Write(LoggingLevel::Default, ss_end.str());
    ShowInputError();
  }
  else
  {
    ss << "The required options are:" << std::endl;
    for (const auto &[arg, options] : mRequiredArguments)
    {
      (void)arg;
      ss << "--" << options.argument << "=\t" << options.description
         << std::endl;
    }
    ss << "The optional arguments are:" << std::endl;
    for (const auto &[arg, options] : mOptionalArguments)
    {
      (void)arg;
      ss << "--" << options.argument << "=\t" << options.description
         << std::endl;
    }
    Logger::Write(LoggingLevel::Default, ss.str());
  }
}

std::string parser::get_value_as_string(const std::string &argument) const
{
  auto throwError = [](const std::string &arg)
  { throw parserException("Option " + arg + " is not defined."); };
  auto arg = to_lower(argument);
  if (auto posReq = mRequiredArguments.find(arg);
      posReq != mRequiredArguments.end())
  {
    auto &optVal = posReq->second.value;
    if (not optVal.has_value())
    {
      throwError(argument);
    }
    return optVal.value();
  }
  else if (auto posOpt = mOptionalArguments.find(arg);
           posOpt != mOptionalArguments.end())
  {
    auto &optVal = posOpt->second.value;
    if (not optVal.has_value())
    {
      throwError(argument);
    }
    return optVal.value();
  }
  else
  {
    throwError(argument);
  }
  return std::string();
}

void parser::add_json_input(const std::string &filename)
{
#ifdef nlohmann_json_FOUND
  using json = nlohmann::json;
  std::ifstream f(filename);
  if (not f.good())
  {
    throw parserException("Can not open the json file " + filename);
  }
  json data = json::parse(f);
  std::vector<KeyValue> input;
  for (const auto &[key, value] : data.items())
  {
    if (StringEndsWith(key, "_comment"))
    {
      continue;
    }
    input.push_back({key, value.get<std::string>()});
  }
  add_input(input);
#else
  throw parserException("nlohmann_json is required.");
#endif
}
void parser::set_help_header(const std::string &header)
{
  mHeader = header;
}

void parser::add_input(const std::vector<KeyValue> &input)
{
  auto throwError = [](const std::string &argument) {
    throw parserException("Argument " + argument + " found but not expected.");
  };
  for (const auto &arg : input)
  {
    auto argument = to_lower(arg.key);
    if (argument == "help")
    {
      print_help();
    }
    else
    {
      if (auto pos = mRequiredArguments.find(argument);
          pos != mRequiredArguments.end())
      {
        pos->second.value = arg.value;
      }
      else if (auto posOptional = mOptionalArguments.find(argument);
               posOptional != mOptionalArguments.end())
      {
        posOptional->second.value = arg.value;
      }
      else
      {
        throwError(argument);
      }
    }
  }
  SetLogger(*this);
}

parser::KeyValue parser::get_key_value(const std::string &input)
{
  auto beginning = input.find_first_not_of("-");
  auto seperator = input.find("=");
  if (seperator == std::string::npos)
  {
    seperator = input.find(" ");
  }
  std::string key   = input.substr(beginning, seperator - beginning);
  std::string value = input.substr(seperator + 1, input.size());
  return {key, value};
}

std::string parser::to_lower(const std::string &input) const
{
  auto arg = input;
  std::transform(arg.begin(),
                 arg.end(),
                 arg.begin(),
                 [](unsigned char c) { return std::tolower(c); });
  return arg;
}

} // namespace BSMPT
