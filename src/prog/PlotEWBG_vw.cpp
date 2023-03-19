// Copyright (C) 2020  Philipp Basler, Margarete Mühlleitner and Jonas
// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

/**
 * @file
 * This program calculates the EWBG eta as a function of vw and varies vw over a
 * given array.
 */

#include <BSMPT/baryo_calculation/CalculateEtaInterface.h>
#include <BSMPT/minimizer/Minimizer.h>
#include <BSMPT/models/ClassPotentialOrigin.h> // for Class_Pot...
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/utility/Logger.h>
#include <BSMPT/utility/parser.h>
#include <BSMPT/utility/utility.h>
#include <algorithm> // for copy, max
#include <fstream>
#include <iostream>
#include <memory>   // for shared_ptr
#include <stdlib.h> // for atof, EXI...
#include <string>   // for operator<<
#include <vector>   // for vector
using namespace std;
using namespace BSMPT;

struct CLIOptions
{
  BSMPT::ModelID::ModelIDs Model{ModelID::ModelIDs::NotSet};
  int Line{};
  std::string InputFile, OutputFile, ConfigFile;
  bool TerminalOutput{false};
  double vw_min{}, vw_max{}, vw_Stepsize{};
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

  const auto SMConstants = GetSMConstants();
  auto argparser         = prepare_parser();
  argparser.add_input(convert_input(argc, argv));
  const CLIOptions args(argparser);
  if (not args.good())
  {
    return EXIT_FAILURE;
  }

  // Set up of BSMPT/Baryo Classes
  Baryo::CalculateEtaInterface EtaInterface(args.ConfigFile, SMConstants);
  std::shared_ptr<Class_Potential_Origin> modelPointer =
      ModelID::FChoose(args.Model, SMConstants);

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
  bool found      = false;

  while (true)
  {
    if (infile.eof()) break;
    std::getline(infile, linestr);
    if (linecounter == 1)
    {
      modelPointer->setUseIndexCol(linestr);
      outfile << linestr << sep;
      outfile << "T_c_var" << sep << "omega_c_var" << sep << "vw_var" << sep
              << "LW_var";
      for (const auto &x : EtaInterface.legend())
        outfile << sep << x + "_var";
      outfile << std::endl;
    }
    else if (linecounter == args.Line)
    {
      modelPointer->initModel(linestr);
      modelPointer->FindSignSymmetries();
      found = true;
      break;
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

  if (args.TerminalOutput) modelPointer->write();
  // CALL: BSMPT-->Phasetransition
  if (args.TerminalOutput)
    Logger::Write(LoggingLevel::ProgDetailed, "PTFinder called...");
  auto EWPT = Minimizer::PTFinder_gen_all(
      modelPointer, 0, 300, args.WhichMinimizer, args.UseMultithreading);

  // SFOEWPT FOUND
  if (EWPT.StatusFlag == Minimizer::MinimizerStatus::SUCCESS and
      C_PT * EWPT.Tc < EWPT.vc)
  {
    if (args.TerminalOutput)
      Logger::Write(LoggingLevel::ProgDetailed, "SFOEWPT found...");
    std::vector<double> vcritical, vbarrier;
    vcritical = EWPT.EWMinimum;
    double TC = EWPT.Tc;
    double vc = EWPT.vc;
    std::vector<double> MinimumPlane;
    // Find the minimum in the symmetric phase. For this minimise at T = Tc + 1
    std::vector<double> vevsymmetricSolution, checksym, startpoint;
    for (std::size_t i = 0; i < modelPointer->get_nVEV(); i++)
      startpoint.push_back(0.5 * vcritical.at(i));
    vevsymmetricSolution = Minimizer::Minimize_gen_all(modelPointer,
                                                       TC + 1,
                                                       checksym,
                                                       startpoint,
                                                       args.WhichMinimizer,
                                                       args.UseMultithreading);
    double vw            = 0;
    if (args.TerminalOutput)
      Logger::Write(LoggingLevel::ProgDetailed,
                    "Currently calculating vw:",
                    __FILE__,
                    __LINE__);
    for (vw = args.vw_min; vw <= args.vw_max; vw += args.vw_Stepsize)
    {
      Logger::Write(LoggingLevel::Default, "\rvw = " + std::to_string(vw));
      auto eta = EtaInterface.CalcEta(vw,
                                      vcritical,
                                      vevsymmetricSolution,
                                      TC,
                                      modelPointer,
                                      args.WhichMinimizer);
      outfile << linestr << sep;
      outfile << TC << sep << vc << sep << vw << sep << EtaInterface.getLW();
      for (auto x : eta)
        outfile << sep << x;
      outfile << std::endl;
    } // END: vw loop
  }   // END: SFOEWPT FOUND
  else
  {
    outfile << -1 << -1 << -1 << -1 << -1 << std::endl;
  } // NO SFOEWPT
  outfile.close();

  return EXIT_SUCCESS;
} // END: Try
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
  Model       = BSMPT::ModelID::getModel(argparser.get_value("model"));
  InputFile   = argparser.get_value("input");
  OutputFile  = argparser.get_value("output");
  Line        = argparser.get_value<int>("line");
  ConfigFile  = argparser.get_value("config");
  vw_min      = argparser.get_value<double>("vw_min");
  vw_max      = argparser.get_value<double>("vw_max");
  vw_Stepsize = argparser.get_value<double>("vw_stepsize");
  try
  {
    TerminalOutput = (argparser.get_value("terminalOutput") == "y");
  }
  catch (BSMPT::parserException &)
  {
    TerminalOutput = false;
  }

  try
  {
    UseGSL = argparser.get_value<bool>("useGSL");
  }
  catch (BSMPT::parserException &)
  {
  }

  try
  {
    UseCMAES = argparser.get_value<bool>("useCMAES");
  }
  catch (BSMPT::parserException &)
  {
  }

  try
  {
    UseNLopt = argparser.get_value<bool>("useNLopt");
  }
  catch (BSMPT::parserException &)
  {
  }

  try
  {
    UseMultithreading = argparser.get_value<bool>("useMultithreading");
  }
  catch (BSMPT::parserException &)
  {
  }

  WhichMinimizer = Minimizer::CalcWhichMinimizer(UseGSL, UseCMAES, UseNLopt);
}

bool CLIOptions::good() const
{
  if (UseGSL and not Minimizer::UseGSLDefault)
  {
    throw std::runtime_error(
        "You set --useGSL=true but GSL was not found during compilation.");
  }
  if (UseCMAES and not Minimizer::UseLibCMAESDefault)
  {
    throw std::runtime_error(
        "You set --useCMAES=true but CMAES was not found during compilation.");
  }
  if (UseNLopt and not Minimizer::UseNLoptDefault)
  {
    throw std::runtime_error(
        "You set --useNLopt=true but NLopt was not found during compilation.");
  }
  if (WhichMinimizer == 0)
  {
    throw std::runtime_error(
        "You disabled all minimizers. You need at least one.");
  }
  if (vw_min <= 0 or vw_min > 1 or vw_max <= 0 or vw_max > 1)
  {
    throw std::runtime_error("The wall velocity has to be between 0 and 1.");
  }
  if (vw_Stepsize == 0)
  {
    throw std::runtime_error("The stepsize has to be larger than 0.");
  }
  if (vw_min > vw_max)
  {
    throw std::runtime_error(
        "The minimal wall velocity has to be smaller than the maximal.");
  }
  if (Model == ModelID::ModelIDs::NotSet)
  {
    Logger::Write(
        LoggingLevel::Default,
        "Your Model parameter does not match with the implemented Models.");
    ShowInputError();
    return false;
  }
  if (Line < 1)
  {
    Logger::Write(LoggingLevel::Default, "Start line counting with 1");
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
  argparser.add_argument(
      "line",
      "The line with the given parameter point. Expects line 1 to "
      "be a legend.",
      true);
  argparser.add_argument("config", "The EWBG config file.", true);

  argparser.add_argument("vw_min", "The minimum wall velocity.", true);
  argparser.add_argument("vw_max", "The maximum wall velocity.", true);
  argparser.add_argument(
      "vw_stepsize", "The stepsize to increase the wall velocity.", true);

  argparser.add_argument(
      "terminalOutput",
      "y/n Turns on additional information in the terminal during "
      "the calculation.",
      false);

  std::stringstream ss;
  ss << "PlotEWBG_vw calculates the EWBG for varying wall velocity "
        "for a given parameter point."
     << std::endl
     << "It is called either by " << std::endl
     << "./PlotEWBG_vw Model Inputfile Outputfile Line vwMin vwStepsize vwMax "
        "EWBGConfigFile"
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
      arguments.emplace_back("--vw_min=" + std::string(argv[5]));
    }
    if (argc >= 7)
    {
      arguments.emplace_back("--vw_stepsize=" + std::string(argv[6]));
    }
    if (argc >= 8)
    {
      arguments.emplace_back("--vw_max=" + std::string(argv[7]));
    }
    if (argc >= 9)
    {
      arguments.emplace_back("--config=" + std::string(argv[8]));
    }
  }
  return arguments;
}
