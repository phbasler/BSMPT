#include <BSMPT/config.h>
#include <BSMPT/utility/Logger.h>
#include <BSMPT/utility/parser.h>
#include <BSMPT/utility/utility.h>
#include <algorithm>
#include <sstream>

#ifdef nlohmann_json_FOUND
#include <fstream>
#include <nlohmann/json.hpp>
#endif

namespace
{

std::string to_lower(const std::string &input)
{
  auto arg = input;
  std::transform(arg.begin(),
                 arg.end(),
                 arg.begin(),
                 [](unsigned char c) { return std::tolower(c); });
  return arg;
}
} // namespace
namespace BSMPT
{

parserException::parserException(const std::string &msg) : message{msg}
{
}

const char *parserException::what() const noexcept
{
  return message.c_str();
}

parser::parser()
{
  add_argument("logginglevel::disabled", "Disable the Logger.", false);
  add_argument("logginglevel::debug",
               "Turn on/off (true/false) the debug output for the logger.",
               false);
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
  add_argument("logginglevel::default",
               "Turn on/off (true/false) the default output.",
               false);

  add_argument("json",
               "Use a json file to define the input instead of additional cli "
               "parameters.",
               false);
}

void parser::enable_minimizer_options()
{
  add_argument("UseGSL",
               "Use the GSL library to minimize the effective potential.",
               false);
  add_argument("UseCMAES",
               "Use the CMAES library to minimize the effective potential",
               false);
  add_argument("UseNLopt",
               "Use the NLopt library to minimize the effective potential",
               false);
  add_argument("UseMultithreading",
               "Enables/Disables multi threading for the minimizers",
               false);
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
  if (required)
  {
    mRequiredArguments.emplace(arg, options);
  }
  else
  {
    mOptionalArguments.emplace(arg, options);
  }
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

std::string parser::get_value(const std::string &argument) const
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

std::string parser::get_value_lower_case(const std::string &argument) const
{
  try
  {
    auto value = get_value(argument);
    return to_lower(value);
  }
  catch (BSMPT::parserException &e)
  {
    throw e;
  }
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
  std::stringstream ss;
  ss << mHeader << std::endl;
  ss << "The required options are:" << std::endl;
  for (const auto &[arg, options] : mRequiredArguments)
  {
    (void)arg;
    ss << "--" << options.argument << ":\t" << options.description << std::endl;
  }
  ss << "The optional arguments are:" << std::endl;
  for (const auto &[arg, options] : mOptionalArguments)
  {
    (void)arg;
    ss << "--" << options.argument << ":\t" << options.description << std::endl;
  }
  Logger::Write(LoggingLevel::Default, ss.str());
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

} // namespace BSMPT
