// Copyright (C) 2020  Philipp Basler, Margarete Mühlleitner and Jonas Müller
// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef PARSER_H
#define PARSER_H

#include <BSMPT/models/IncludeAllModels.h>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

namespace BSMPT
{

/**
 * @brief ShowInputError shows all the available models in the terminal
 */
void ShowInputError();

class parserException : public std::exception
{
public:
  virtual const char *what() const noexcept;
  parserException(const std::string &msg);

private:
  std::string message;
};

/**
 * @brief The parser class provides the argument parser for the CLI and JSON
 * methods. This is case insensitive.
 */
class parser
{
public:
  parser();
  parser(const bool &enable_column_output);
  /**
   * @brief enable_minimizer_options Enables the options regarding the
   * minimizer.
   */
  void enable_minimizer_options();
  /**
   * @brief add_argument Silently adds a new argument for the parser options.
   * @param argument The parameter name for the CLI or JSON key. This will be
   * treated case insensitive.
   * @param description The description for the parameter shown in help.
   * @param default_val The default value for optional parameters.
   */
  void add_argument_only_display(const std::string &argument,
                                 const std::string &description,
                                 const std::string &default_val);
  /**
   * @brief add_argument Silently adds a new argument for the parser options.
   * @param argument The parameter name for the CLI or JSON key. This will be
   * treated case insensitive.
   * @param required Decide if it is a required parameter or not.
   */
  void add_argument(const std::string &argument, bool required);
  /**
   * @brief add_argument Adds a new argument for the parser options.
   * @param argument The parameter name for the CLI or JSON key. This will be
   * treated case insensitive.
   * @param description The description for the parameter shown in help.
   * @param required Decide if it is a required parameter or not.
   */
  void add_argument(const std::string &argument,
                    const std::string &description,
                    bool required);
  /**
   * @brief add_argument Adds a new argument for the parser options.
   * @param argument The parameter name for the CLI or JSON key. This will be
   * treated case insensitive.
   * @param description The description for the parameter shown in help.
   * @param default_val The default value for optional parameters.
   * @param required Decide if it is a required parameter or not.
   */
  void add_argument(const std::string &argument,
                    const std::string &description,
                    const std::string &default_val,
                    bool required);
  /**
   * @brief add_subtext add subtext to description column
   * @param subtext string to print in description column
   */
  void add_subtext(const std::string &subtext);
  /**
   * @brief add_input Add the vector with each input as it is given in the CLI
   * in the form "--argument=value".
   * @param input The vector containing the CLI inputs.
   * @throws parserException if the argument was not expected by the parser.
   */
  void add_input(const std::vector<std::string> &input);
  /**
   * @brief print_help Prints the header and the arguments with their
   * description.
   */
  void print_help() const;
  /**
   * @brief get_value Get the value for the required parameter.
   * @param argument The required CLI name or JSON key.
   * @return The value if the parameter was set.
   * @throws parserException if the argument was not set or if the value can not
   * be casted from a string to the type.
   */
  template <typename T = std::string>
  T get_value(const std::string &argument) const
  {
    auto value = get_value_as_string(argument);
    if constexpr (std::is_same<T, std::string>::value)
    {
      return value;
    }
    else if constexpr (std::is_same<T, int>::value)
    {
      return std::stoi(value);
    }
    else if constexpr (std::is_same<T, unsigned int>::value)
    {
      return std::stoul(value);
    }
    else if constexpr (std::is_same<T, double>::value)
    {
      return std::stod(value);
    }
    else if constexpr (std::is_same<T, bool>::value)
    {
      auto lower = to_lower(value);
      return lower == "true";
    }
    else
    {
      return static_cast<T>(value);
    }
  }

  /**
   * @brief all_required_set Check ifs all required parameter are set.
   * @return true/false if all required parameters are set.
   */
  bool all_required_set() const;
  /**
   * @brief check_required_parameters Calls all_required_set() and throws a
   * parserException if not all required parameters are set.
   * @throws parserException if not all required parameter are set.
   */
  void check_required_parameters() const;
  /**
   * @brief set_help_header Sets the header which is printed in print_help()
   * @param header The header to print.
   */
  void set_help_header(const std::string &header);

private:
  struct Options
  {
    std::string argument;
    std::string description;
    std::optional<std::string> value;
    std::string default_val = "";
  };
  struct KeyValue
  {
    std::string key;
    std::string value;
    KeyValue() = default;
    KeyValue(const std::string &Key, const std::string &Value)
        : key{Key}
        , value{Value}
    {
    }
  };

  std::vector<std::pair<std::string, Options>> mOrderedArguments;
  std::unordered_map<std::string, Options> mRequiredArguments;
  std::unordered_map<std::string, Options> mOptionalArguments;

  bool extra_column_output = false;

  /**
   * @brief get_value Get the value for the required parameter.
   * @param argument The required CLI name or JSON key.
   * @return The value if the parameter was set.
   * @throws parserException if the argument was not set or if the value can not
   * be casted from a string to the type.
   */
  std::string get_value_as_string(const std::string &argument) const;
  /**
   * @brief add_json_input Add the input parameters with a json file.
   * @param filename The file containing the json.
   */
  void add_json_input(const std::string &filename);
  /**
   * @brief add_input Set the key value pairs.
   * @param input The vector containing all key value pairs.
   */
  void add_input(const std::vector<KeyValue> &input);
  /**
   * @brief get_key_value Converts the CLI input in the form "--argName=value"
   * into a key value pair.
   * @param input The CLI option in the form of "--argName=foo".
   * @return Returns the matching KeyValue.
   */

  std::string to_lower(const std::string &input) const;

  KeyValue get_key_value(const std::string &input);
  std::string mHeader;

  /**
   * @brief Helper to avoid multiple prints of the help output
   */
  mutable bool mHelpAlreadyPrinted{false};
};
} // namespace BSMPT
#endif // PARSER_H
