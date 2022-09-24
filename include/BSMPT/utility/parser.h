#ifndef PARSER_H
#define PARSER_H

#include <string>
#include <unordered_map>

namespace BSMPT
{

class parserException : public std::exception
{
public:
  char *what();
  parserException(const std::string &msg);

private:
  std::string message;
};

/**
 * @brief The parser class provides the base class for parsing CLI options
 */
class parser
{
public:
  parser();
  void add_argument(const std::string &argument,
                    const std::string &description);
  void add_input(int argc, char *argv[]);
  void print_help() const;
  std::string get_value(const std::string &argument) const;

private:
  struct Options
  {
    std::string argument;
    std::string description;
    std::string value;
  };
  std::unordered_map<std::string, Options> mArguments;
  void add_cli_input(int argc, char *argv[]);
  void add_json_input(const std::string &filename);
};
} // namespace BSMPT
#endif // PARSER_H
