// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include <BSMPT/utility/Logger.h>
#include <BSMPT/utility/parser.h>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <sstream>

using Approx = Catch::Approx;

TEST_CASE("Ask unexpected option", "[parser]")
{
  auto parser = BSMPT::parser();
  REQUIRE_THROWS_AS(parser.get_value("non existent"), BSMPT::parserException);
}

TEST_CASE("Check bool input", "[parser]")
{
  auto parser = BSMPT::parser();
  std::string argName{"arg"}, argNameF{"arg2"};
  parser.add_argument(argName, "foo", true);
  parser.add_argument(argNameF, "foo", true);
  std::vector<std::string> input;
  input.emplace_back("--" + argName + "=True");
  input.emplace_back("--" + argNameF + "=foo");
  parser.add_input(input);
  REQUIRE(parser.get_value<bool>(argName));
  REQUIRE_FALSE(parser.get_value<bool>(argNameF));
}

TEST_CASE("Check int input", "[parser]")
{
  auto parser = BSMPT::parser();
  int value   = 3;
  std::string argName{"arg"};
  parser.add_argument(argName, "foo", true);
  std::vector<std::string> input;
  input.emplace_back("--" + argName + "=" + std::to_string(value));
  parser.add_input(input);
  REQUIRE(parser.get_value<int>(argName) == value);
}

TEST_CASE("Check double input", "[parser]")
{
  auto parser  = BSMPT::parser();
  double value = 3.2;
  std::string argName{"arg"};
  parser.add_argument(argName, "foo", true);
  std::vector<std::string> input;
  input.emplace_back("--" + argName + "=" + std::to_string(value));
  parser.add_input(input);
  REQUIRE(parser.get_value<double>(argName) == value);
}

TEST_CASE("Set unexpected input", "[parser]")
{
  auto parser = BSMPT::parser();
  std::string argName{"ExampleArg"};
  std::vector<std::string> input;
  input.emplace_back("--" + argName + "=foo");
  REQUIRE_THROWS_AS(parser.add_input(input), BSMPT::parserException);
}

TEST_CASE("Add argument and retrieve it", "[parser]")
{
  auto parser = BSMPT::parser();
  std::string argName{"ExampleArg"};
  parser.add_argument(argName, "some example description", false);
  std::vector<std::string> input;
  input.emplace_back("--" + argName + "=foo");
  parser.add_input(input);
  REQUIRE_NOTHROW(parser.get_value(argName));
}

TEST_CASE("Check for logger print", "[parser]")
{
  auto parser = BSMPT::parser();
  std::string argName{"ExampleArg"};
  std::string description{"Example description"};
  std::string header{"Some header"};
  parser.add_argument(argName, description, false);
  parser.set_help_header(header);
  std::stringstream ss;
  BSMPT::Logger::SetOStream(ss);
  parser.print_help();
  auto str = ss.str();
  CHECK(str.find(header) != std::string::npos);
  CHECK(str.find(argName) != std::string::npos);
  CHECK(str.find(description) != std::string::npos);
  BSMPT::Logger::SetOStream(std::cout);
}

TEST_CASE("Check for all required parameters", "[parser]")
{
  auto parser = BSMPT::parser();
  std::string argName{"ExampleArg"};
  std::string description{"Example description"};
  std::string header{"Some header"};
  parser.add_argument(argName, description, true);
  REQUIRE_THROWS_AS(parser.check_required_parameters(), BSMPT::parserException);
}

TEST_CASE("Check parser input error is disabled on Loglevel None", "[parser]")
{
  using namespace BSMPT;
  std::stringstream ss;
  Logger::SetOStream(ss);
  Logger::Disable();
  ShowInputError();
  Logger::Write(LoggingLevel::Default, "Some output");
  std::string output = ss.str();
  REQUIRE(output.empty());
  Logger::SetOStream(std::cout);
  Logger::RestoreDefaultLevels();
}