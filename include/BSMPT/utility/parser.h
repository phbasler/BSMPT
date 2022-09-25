#ifndef PARSER_H
#define PARSER_H

#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

namespace BSMPT
{

class parserException : public std::exception
{
public:
  virtual const char *what() const noexcept;
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
  void enable_minimizer_options();
  void add_argument(const std::string &argument,
                    const std::string &description,
                    bool required);
  void add_input(const std::vector<std::string> &input);
  void print_help() const;
  std::string get_value(const std::string &argument) const;
  std::string get_value_lower_case(const std::string &argument) const;
  bool all_required_set() const;
  void check_required_parameters() const;
  void set_help_header(const std::string &header);

private:
  struct Options
  {
    std::string argument;
    std::string description;
    std::optional<std::string> value;
  };
  std::unordered_map<std::string, Options> mRequiredArguments;
  std::unordered_map<std::string, Options> mOptionalArguments;
  void add_cli_input(const std::vector<std::string> &input);
  void add_json_input(const std::string &filename);
  std::string mHeader;
};
} // namespace BSMPT
#endif // PARSER_H
