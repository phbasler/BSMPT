// Copyright (C) 2020  Philipp Basler, Margarete Mühlleitner and Jonas
// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

/**
 * @file
 * Calculates the Counterterms for mulitple points and adds them at the end of
 * the line. The sequence is the same as used in Potential::set_CT_Pot_Par
 */

#include <BSMPT/models/ClassPotentialOrigin.h> // for Class_Potential_Origin
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/utility/Logger.h>
#include <BSMPT/utility/parser.h>
#include <BSMPT/utility/utility.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>   // for unique_ptr
#include <stdlib.h> // for atoi, EXIT_FAILURE
#include <string>   // for operator<<, getline
#include <utility>  // for pair
#include <vector>   // for vector
using namespace std;
using namespace BSMPT;

struct CLIOptions
{
public:
  BSMPT::ModelID::ModelIDs Model{ModelID::ModelIDs::NotSet};
  int FirstLine{}, LastLine{};
  std::string InputFile, OutputFile;
  bool TerminalOutput{false};

  CLIOptions(const BSMPT::parser &argparser);
  bool good() const;

private:
  BSMPT::parser mArgparser;
};

BSMPT::parser prepare_parser();

std::vector<std::string> convert_input(int argc, char *argv[]);

int main(int argc, char *argv[])
try
{

  auto argparser = prepare_parser();
  argparser.add_input(convert_input(argc, argv));
  const CLIOptions args(argparser);

  if (not args.good())
  {
    return EXIT_FAILURE;
  }

  int linecounter = 1;
  std::ifstream infile(args.InputFile);
  if (!infile.good())
  {
    Logger::Write(LoggingLevel::Default, "Input file not found ");
    return EXIT_FAILURE;
  }
  std::ofstream outfile(args.OutputFile);
  if (!outfile.good())
  {
    Logger::Write(LoggingLevel::Default,
                  "Can not create file " + args.OutputFile);
    return EXIT_FAILURE;
  }
  std::string linestr;

  std::unique_ptr<Class_Potential_Origin> modelPointer =
      ModelID::FChoose(args.Model);

  while (getline(infile, linestr))
  {
    if (linecounter > args.LastLine) break;
    if (linecounter == 1)
    {
      modelPointer->setUseIndexCol(linestr);
      outfile << linestr;
      outfile << sep << modelPointer->addLegendCT();
      outfile << std::endl;
    }

    if (linecounter >= args.FirstLine and linecounter <= args.LastLine and
        linecounter != 1)
    {
      std::pair<std::vector<double>, std::vector<double>> parameters =
          modelPointer->initModel(linestr);
      outfile << linestr << sep << parameters.second << std::endl;
    }
    linecounter++;
    if (infile.eof()) break;
  }
  outfile.close();
  return EXIT_SUCCESS;
}
catch (int)
{
  return EXIT_SUCCESS;
}
catch (exception &e)
{
  Logger::Write(LoggingLevel::Default, e.what());
  return EXIT_FAILURE;
}

CLIOptions::CLIOptions(const BSMPT::parser &argparser)
{
  argparser.check_required_parameters();
  Model      = BSMPT::ModelID::getModel(argparser.get_value("model"));
  InputFile  = argparser.get_value("input");
  OutputFile = argparser.get_value("output");
  FirstLine  = argparser.get_value<int>("FirstLine");
  LastLine   = argparser.get_value<int>("LastLine");
  try
  {
    TerminalOutput = (argparser.get_value("TerminalOutput") == "y");
  }
  catch (BSMPT::parserException &)
  {
    TerminalOutput = false;
  }
}

bool CLIOptions::good() const
{

  if (Model == ModelID::ModelIDs::NotSet)
  {
    Logger::Write(
        LoggingLevel::Default,
        "Your Model parameter does not match with the implemented Models.");
    ShowInputError();
    return false;
  }
  if (FirstLine < 1)
  {
    Logger::Write(LoggingLevel::Default, "Start line counting with 1");
    return false;
  }
  if (FirstLine > LastLine)
  {
    Logger::Write(LoggingLevel::Default, "Firstline is smaller then LastLine ");
    return false;
  }
  return true;
}

BSMPT::parser prepare_parser()
{
  BSMPT::parser argparser;
  argparser.add_argument("model", "The model you want to investigate.", true);
  argparser.add_argument("input", "The input file in tsv format.", true);
  argparser.add_argument("output", "The output file in tsv format.", true);
  argparser.add_argument("FirstLine",
                         "The first line in the input file to calculate the "
                         "EWPT. Expects line 1 to be a legend.",
                         true);
  argparser.add_argument(
      "LastLine",
      "The last line in the input file to calculate the EWPT.",
      true);
  argparser.add_argument(
      "TerminalOutput",
      "y/n Turns on additional information in the terminal during "
      "the calculation.",
      false);

  std::stringstream ss;
  ss << "CalcCT calculates the counterterm parameters for the given "
        "paramter points"
     << std::endl
     << "It is called either by " << std::endl
     << "./CalcCT model input output FirstLine LastLine" << std::endl
     << "or with the following arguments" << std::endl;
  argparser.set_help_header(ss.str());

  return argparser;
}

std::vector<std::string> convert_input(int argc, char *argv[])
{
  std::vector<std::string> arguments;
  if (argc == 1) return arguments;
  auto first_arg = std::string(argv[1]);

  bool UsePrefix =
      StringStartsWith(first_arg, "--") or StringStartsWith(first_arg, "-");

  if (UsePrefix)
  {
    for (int i{1}; i < argc; ++i)
    {
      arguments.emplace_back(argv[i]);
    }
  }
  else
  {
    if (argc >= 2)
    {
      arguments.emplace_back("--model=" + std::string(argv[1]));
    }
    if (argc >= 3)
    {
      arguments.emplace_back("--input=" + std::string(argv[2]));
    }
    if (argc >= 4)
    {
      arguments.emplace_back("--output=" + std::string(argv[3]));
    }
    if (argc >= 5)
    {
      arguments.emplace_back("--FirstLine=" + std::string(argv[4]));
    }
    if (argc >= 6)
    {
      arguments.emplace_back("--LastLine=" + std::string(argv[5]));
    }
  }
  return arguments;
}
