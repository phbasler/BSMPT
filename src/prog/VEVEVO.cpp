// Copyright (C) 2020  Philipp Basler, Margarete Mühlleitner and Jonas
// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

/**
 * @file
 * This program calculates the development of the VeVs with the Temperature.
 */

#include <BSMPT/minimizer/Minimizer.h>
#include <BSMPT/models/ClassPotentialOrigin.h> // for Class_Potential_Origin
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/utility/Logger.h>
#include <BSMPT/utility/parser.h>
#include <BSMPT/utility/utility.h>
#include <algorithm> // for copy, max
#include <fstream>
#include <iomanip> // for operator<<, setprecision
#include <iostream>
#include <math.h>   // for sqrt, abs
#include <memory>   // for shared_ptr, __shared_...
#include <stdlib.h> // for atof, EXIT_FAILURE, atoi
#include <string>   // for getline, operator<<
#include <utility>  // for pair
#include <vector>   // for vector
using namespace std;
using namespace BSMPT;

struct CLIOptions
{
  BSMPT::ModelID::ModelIDs Model{ModelID::ModelIDs::NotSet};
  int Line{};
  std::string InputFile, OutputFile;
  double TemperatureStart{}, TemperatureStep{}, TemperatureEnd{};
  bool UseGSL{Minimizer::UseGSLDefault};
  bool UseCMAES{Minimizer::UseLibCMAESDefault};
  bool UseNLopt{Minimizer::UseNLoptDefault};
  int WhichMinimizer{Minimizer::WhichMinimizerDefault};
  bool UseMultithreading{true};

  CLIOptions(const BSMPT::parser &argparser);
  bool good() const;
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

  std::vector<double> sol, start, solPot;

  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(args.Model);
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
  int linecounter = 1;

  bool found = false;
  while (true)
  {
    if (infile.eof()) break;
    std::getline(infile, linestr);
    if (linecounter == 1)
    {
      modelPointer->setUseIndexCol(linestr);
    }
    else if (linecounter == args.Line)
    {
      std::pair<std::vector<double>, std::vector<double>> parameters =
          modelPointer->initModel(linestr);
      found = true;
    }
    else if (linecounter > args.Line)
      break;
    linecounter++;
    if (infile.eof()) break;
  }
  infile.close();
  if (!found)
  {
    Logger::Write(LoggingLevel::Default, "Line not found !");
    return EXIT_FAILURE;
  }

  std::vector<double> Check;
  double vev{0.0};

  outfile << std::setprecision(16);

  outfile << "T" << sep << "v";
  for (auto x : modelPointer->addLegendVEV())
    outfile << sep << x;
  outfile << sep << "Veff(v,T)" << std::endl;

  for (double Temp = args.TemperatureStart; Temp <= args.TemperatureEnd;
       Temp += args.TemperatureStep)
  {
    start.clear();
    if (Temp == args.TemperatureStart)
    {
      start = modelPointer->get_vevTreeMin();
    }
    else
    {
      start = sol;
    }
    sol.clear();
    Check.clear();
    solPot.clear();
    sol    = Minimizer::Minimize_gen_all(modelPointer,
                                      Temp,
                                      Check,
                                      start,
                                      args.WhichMinimizer,
                                      args.UseMultithreading);
    solPot = modelPointer->MinimizeOrderVEV(sol);
    vev    = modelPointer->EWSBVEV(solPot);

    outfile << Temp << sep;
    outfile << vev << sep;
    outfile << sol;
    outfile << sep << modelPointer->VEff(solPot, Temp, 0);
    outfile << std::endl;
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

bool CLIOptions::good() const
{
  if (UseGSL and not Minimizer::UseGSLDefault)
  {
    throw std::runtime_error(
        "You set --UseGSL=true but GSL was not found during compilation.");
  }
  if (UseCMAES and not Minimizer::UseLibCMAESDefault)
  {
    throw std::runtime_error(
        "You set --UseCMAES=true but CMAES was not found during compilation.");
  }
  if (UseNLopt and not Minimizer::UseNLoptDefault)
  {
    throw std::runtime_error(
        "You set --UseNLopt=true but NLopt was not found during compilation.");
  }
  if (WhichMinimizer == 0)
  {
    throw std::runtime_error(
        "You disabled all minimizers. You need at least one.");
  }
  if (TemperatureStep <= 0)
  {
    throw std::runtime_error("The stepsize has to be larger than 0.");
  }
  if (TemperatureStart < 0)
  {
    throw std::runtime_error(
        "The starting value of the temperature can not be negative.");
  }
  if (TemperatureEnd < TemperatureStart)
  {
    throw std::runtime_error("The minimal value for the temperature is lower "
                             "then the maximal value.");
  }
  if (Model == ModelID::ModelIDs::NotSet)
  {

    Logger::Write(
        LoggingLevel::Default,
        "Your Model parameter does not match with the implemented Models.");
    ShowInputError();
    return false;
  }
  return true;
}

CLIOptions::CLIOptions(const BSMPT::parser &argparser)
{
  argparser.check_required_parameters();
  Model            = BSMPT::ModelID::getModel(argparser.get_value("model"));
  InputFile        = argparser.get_value("input");
  OutputFile       = argparser.get_value("output");
  Line             = std::stoi(argparser.get_value("line"));
  TemperatureStart = std::stod(argparser.get_value("TemperatureStart"));
  TemperatureEnd   = std::stod(argparser.get_value("TemperatureEnd"));
  TemperatureStep  = std::stod(argparser.get_value("TemperatureStep"));

  try
  {
    UseGSL = argparser.get_value_lower_case("UseGSL") == "true";
  }
  catch (BSMPT::parserException &)
  {
  }

  try
  {
    UseCMAES = argparser.get_value_lower_case("UseCMAES") == "true";
  }
  catch (BSMPT::parserException &)
  {
  }

  try
  {
    UseNLopt = argparser.get_value_lower_case("UseNLopt") == "true";
  }
  catch (BSMPT::parserException &)
  {
  }

  try
  {
    UseMultithreading =
        argparser.get_value_lower_case("UseMultithreading") == "true";
  }
  catch (BSMPT::parserException &)
  {
  }

  WhichMinimizer = Minimizer::CalcWhichMinimizer(UseGSL, UseCMAES, UseNLopt);
}

BSMPT::parser prepare_parser()
{
  BSMPT::parser argparser;
  argparser.add_argument("model", "The model you want to investigate.", true);
  argparser.add_argument("input", "The input file in tsv format.", true);
  argparser.add_argument("output", "The output file in tsv format.", true);
  argparser.add_argument(
      "line",
      "The line with the given parameter point. Expects line 1 to "
      "be a legend.",
      true);

  argparser.add_argument(
      "TemperatureStart",
      "The starting temperature to calculate the global minimum.",
      true);
  argparser.add_argument(
      "TemperatureStep", "The stepsize for the temperature.", true);
  argparser.add_argument(
      "TemperatureEnd",
      "The last temperature to calculate the global minimum.",
      true);
  argparser.add_argument(
      "TerminalOutput",
      "y/n Turns on additional information in the terminal during "
      "the calculation.",
      false);

  argparser.add_argument(
      "vw",
      "Wall velocity for the EWBG calculation. Default value of 0.1.",
      false);

  std::stringstream ss;
  ss << "VEVEVO calculates the evolution of the global minimum with "
        "rising temperature for a given parameter point"
     << std::endl
     << "It is called either by " << std::endl
     << "./VEVEVO Model Inputfile Outputfile Line TemperatureStart "
        "TemperatureStep TemperatureEnd"
     << std::endl
     << "or with the following arguments" << std::endl;
  argparser.set_help_header(ss.str());

  argparser.enable_minimizer_options();

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
      arguments.emplace_back("--line=" + std::string(argv[4]));
    }
    if (argc >= 6)
    {
      arguments.emplace_back("--TemperatureStart=" + std::string(argv[5]));
    }
    if (argc >= 7)
    {
      arguments.emplace_back("--TemperatureStep=" + std::string(argv[6]));
    }
    if (argc >= 8)
    {
      arguments.emplace_back("--TemperatureEnd=" + std::string(argv[7]));
    }
  }
  return arguments;
}
