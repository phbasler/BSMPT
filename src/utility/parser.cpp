#include <BSMPT/utility/Logger.h>
#include <BSMPT/utility/parser.h>
#include <sstream>

namespace
{
std::tuple<std::string, std::string> get_key_value(const std::string &input)
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
} // namespace
namespace BSMPT
{

parserException::parserException(const std::string &msg) : message(msg)
{
}

char *parserException::what()
{
  return const_cast<char *>(message.c_str());
}

parser::parser()
{
  add_argument("--logginglevel::disabled", "Disable the Logger.");
  add_argument("--logginglevel::debug",
               "Turn on/off (true/false) the debug output for the logger.");
  add_argument("--logginglevel::ewbgdetailed",
               "Turn on/off (true/false) the output for the EWBG calculation.");
  add_argument(
      "--logginglevel::progdetailed",
      "Turn on/off (true/false) the additional output for the main program.");
  add_argument(
      "--logginglevel::minimizerdetailed",
      "Turn on/off (true/false) the additional output of the minimizers.");
}

void parser::add_argument(const std::string &argument,
                          const std::string &description)
{
  Options options;
  options.argument    = argument;
  options.description = description;
  options.value       = std::string();
  mArguments.emplace(argument, options);
}

void parser::add_input(int argc, char *argv[])
{
  if (argc == 2)
  {
    std::string input(argv[1]);
    auto [argument, value] = get_key_value(input);
    if (argument == "help")
    {
      print_help();
      return;
    }
    else if (argument == "config")
    {
      add_json_input(value);
      return;
    }
  }

  add_cli_input(argc, argv);
}

void parser::add_cli_input(int argc, char *argv[])
{
  for (int i{1}; i < argc; ++i)
  {
    std::string input(argv[i]);
    auto [argument, value] = get_key_value(input);
    if (argument == "help")
    {
      print_help();
    }
    else
    {
      if (auto pos = mArguments.find(argument); pos != mArguments.end())
      {
        pos->second.value = value;
      }
      else
      {
        throw parserException("Argument " + argument +
                              " found but not expected.");
      }
    }
  }
}

std::string parser::get_value(const std::string &argument) const
{
  auto pos = mArguments.find(argument);
  if (pos == mArguments.end())
  {
    throw parserException("Option " + argument + " is not defined.");
  }
  return pos->second.value;
}

void parser::print_help() const
{
  std::stringstream ss;
  ss << "The available options are:" << std::endl;
  for (const auto &[arg, options] : mArguments)
  {
    (void)arg;
    ss << options.argument << ":\t" << options.description << std::endl;
  }
  Logger::Write(LoggingLevel::Default, ss.str());
}

void parser::add_json_input(const std::string &filename)
{
  (void)filename;
}

} // namespace BSMPT
