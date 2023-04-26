// Copyright (C) 2020  Philipp Basler, Margarete Mühlleitner and Jonas
// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

/**
 * @file
 * Calculates the electroweak baryogenesis for a given Inputfile for a given
 * subset of lines in the file and adds it at the end of the line in the format
 * T_c v_c all single vevs. One parameter point per line.
 *
 */

#include <BSMPT/baryo_calculation/CalculateEtaInterface.h>
#include <BSMPT/baryo_calculation/transport_equations.h> // for GSL_integ...
#include <BSMPT/minimizer/Minimizer.h>
#include <BSMPT/models/ClassPotentialOrigin.h> // for Class_Pot...
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/utility/Logger.h>
#include <BSMPT/utility/parser.h>
#include <BSMPT/utility/utility.h>
#include <algorithm> // for max, copy
#include <fstream>
#include <iostream>
#include <memory>   // for shared_ptr
#include <stdlib.h> // for atoi, std::size_t
#include <string>   // for string
#include <vector>   // for vector

using namespace std;
using namespace BSMPT;

struct CLIOptions
{
  BSMPT::ModelID::ModelIDs Model{ModelID::ModelIDs::NotSet};
  int FirstLine{}, LastLine{};
  std::string InputFile, OutputFile, ConfigFile;
  bool TerminalOutput{false};
  double vw{0.1};
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

  auto argparser = prepare_parser();
  argparser.add_input(convert_input(argc, argv));
  const CLIOptions args(argparser);

  if (not args.good())
  {
    return EXIT_FAILURE;
  }

  // Init: Interface Class for the different transport methods
  Baryo::CalculateEtaInterface EtaInterface(args.ConfigFile, SMConstants);

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
  std::shared_ptr<Class_Potential_Origin> modelPointer = ModelID::FChoose(
      args.Model,
      SMConstants); // Declare the model pointer with the necessary parameters
  std::vector<std::string> etaLegend =
      EtaInterface.legend(); // Declare the vector for the PTFinder algorithm
                             // Begin: Input Read
  while (getline(infile, linestr))
  {
    if (linecounter > args.LastLine)
      break;
    else if (linecounter == 1)
    {
      // Write legend
      modelPointer->setUseIndexCol(linestr);
      outfile << linestr;
      for (const auto &x : modelPointer->addLegendCT())
        outfile << sep << x + "_EWBG";
      for (const auto &x : modelPointer->addLegendTemp())
        outfile << sep << x + "_EWBG";
      outfile << sep << "vw";
      outfile << sep << "L_W";
      outfile << sep << "top_sym_phase";
      outfile << sep << "top_brk_phase";
      outfile << sep << "bot_sym_phase";
      outfile << sep << "bot_brk_phase";
      outfile << sep << "tau_sym_phase";
      outfile << sep << "tau_brk_phase";
      outfile << sep << etaLegend;
      outfile << std::endl;
    }
    else if (linecounter >= args.FirstLine and linecounter <= args.LastLine and
             linecounter != 1)
    {
      if (args.TerminalOutput)
      {
        Logger::Write(LoggingLevel::ProgDetailed,
                      "Currently at line " + std::to_string(linecounter));
      }
      // Begin: Parameter Set Up for BSMPT
      auto parameters = modelPointer->initModel(linestr);
      modelPointer->FindSignSymmetries();
      if (args.FirstLine == args.LastLine)
      {
        modelPointer->write();
        Logger::Write(LoggingLevel::Default, "vw = " + std::to_string(args.vw));
      }
      if (args.TerminalOutput)
        Logger::Write(LoggingLevel::ProgDetailed, "Calling PTFinder");

      // Call: BSMPT
      auto EWPT = Minimizer::PTFinder_gen_all(
          modelPointer, 0, 300, args.WhichMinimizer, args.UseMultithreading);
      // Define parameters for eta
      std::vector<double> eta;
      if (EWPT.StatusFlag == Minimizer::MinimizerStatus::SUCCESS and
          C_PT * EWPT.Tc < EWPT.vc)
      {
        if (args.TerminalOutput)
          Logger::Write(LoggingLevel::ProgDetailed, "SFOEWPT found...");
        // Find the minimum in the symmetric phase. For this minimise at T = Tc
        // + 1
        std::vector<double> vevsymmetricSolution, checksym, startpoint;
        for (const auto &el : EWPT.EWMinimum)
          startpoint.push_back(0.5 * el);
        vevsymmetricSolution =
            Minimizer::Minimize_gen_all(modelPointer,
                                        EWPT.Tc + 1,
                                        checksym,
                                        startpoint,
                                        args.WhichMinimizer,
                                        args.UseMultithreading);
        // Call: Calculation of eta in the different implemented approaches
        if (args.TerminalOutput)
          Logger::Write(LoggingLevel::ProgDetailed, "Calling CalcEta...");
        eta = EtaInterface.CalcEta(args.vw,
                                   EWPT.EWMinimum,
                                   vevsymmetricSolution,
                                   EWPT.Tc,
                                   modelPointer,
                                   args.WhichMinimizer);
        // Outfile
        outfile << linestr;
        outfile << sep << parameters.second;
        outfile << sep << EWPT.Tc << sep << EWPT.vc;
        outfile << sep << EWPT.vc / EWPT.Tc;
        outfile << sep << EWPT.EWMinimum;
        outfile << sep << args.vw;
        outfile << sep << EtaInterface.getLW();
        outfile << sep << EtaInterface.getSymmetricCPViolatingPhase_top();
        outfile << sep << EtaInterface.getBrokenCPViolatingPhase_top();
        outfile << sep << EtaInterface.getSymmetricCPViolatingPhase_bot();
        outfile << sep << EtaInterface.getBrokenCPViolatingPhase_bot();
        outfile << sep << EtaInterface.getSymmetricCPViolatingPhase_tau();
        outfile << sep << EtaInterface.getBrokenCPViolatingPhase_tau();
        outfile << sep << eta;
        outfile << std::endl;
      } // END: SFOEWPT found
      else
      { // No SFOEWPT provided
        outfile << linestr;
        outfile << sep << parameters.second;
        outfile << sep << EWPT.Tc << sep << EWPT.vc;
        outfile << sep << EWPT.EWMinimum;
        outfile << sep << EWPT.vc / EWPT.Tc;
        outfile << sep << args.vw;
        outfile << sep << -1;  // LW
        outfile << sep << -50; // top sym CP phase
        outfile << sep << -50; // top brk CP phase
        outfile << sep << -50; // bot sym CP phase
        outfile << sep << -50; // bot brk CP phase
        outfile << sep << -50; // tau sym CP phase
        outfile << sep << -50; // tau brk CP phase
        for (std::size_t i = 0; i < etaLegend.size(); i++)
          outfile << sep << 0;
        outfile << std::endl;
      } // END: No SFOEWPT

      if (args.FirstLine == args.LastLine)
      {
        std::stringstream ss;
        auto dimensionnames = modelPointer->addLegendTemp();
        ss << "Succeded ? " << static_cast<int>(EWPT.StatusFlag) << sep
           << " (1 = Success , -1 = v/T reached a value below " << C_PT
           << " during the calculation) \n";
        if (EWPT.StatusFlag == Minimizer::MinimizerStatus::SUCCESS)
        {
          ss << std::scientific;
          ss << dimensionnames.at(1) << " = " << EWPT.vc << " GeV\n";
          ss << dimensionnames.at(0) << " = " << EWPT.Tc << " GeV\n";
          ss << "xi_c = " << dimensionnames.at(2) << " = " << EWPT.vc / EWPT.Tc
             << std::endl;
          for (std::size_t i = 0; i < modelPointer->get_nVEV(); i++)
          {
            ss << dimensionnames.at(i + 3) << " = " << EWPT.EWMinimum.at(i)
               << " GeV\n";
          }
          ss << "The Wall thickness is given by L_W  = " << EtaInterface.getLW()
             << "GeV^-2\n"
             << "L_W * T = " << EtaInterface.getLW() * EWPT.Tc << "\n";
          for (std::size_t i = 0; i < etaLegend.size(); i++)
            ss << etaLegend.at(i) << " = " << eta.at(i) << std::endl;
        }
        Logger::Write(LoggingLevel::Default, ss.str());

      } // END: LineStart == LineEnd
    }   // END: Valid Line
    linecounter++;
    if (infile.eof()) break;
  } // END: Input Read
  // Closing & Free
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
  Model      = BSMPT::ModelID::getModel(argparser.get_value("model"));
  InputFile  = argparser.get_value("input");
  OutputFile = argparser.get_value("output");
  FirstLine  = argparser.get_value<int>("firstLine");
  LastLine   = argparser.get_value<int>("lastLine");
  ConfigFile = argparser.get_value("config");
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

  try
  {
    vw = std::stod(argparser.get_value("vw"));
  }
  catch (BSMPT::parserException &)
  {
  }
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
  if (vw <= 0 or vw > 1)
  {
    throw std::runtime_error("The wall velocity has to be between 0 and 1.");
  }
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
  argparser.add_argument("firstLine",
                         "The first line in the input file to calculate the "
                         "EWPT. Expects line 1 to be a legend.",
                         true);
  argparser.add_argument(
      "lastLine",
      "The last line in the input file to calculate the EWPT.",
      true);
  argparser.add_argument("config", "The EWBG config file.", true);
  argparser.add_argument(
      "terminalOutput",
      "y/n Turns on additional information in the terminal during "
      "the calculation.",
      false);

  argparser.add_argument(
      "vw",
      "Wall velocity for the EWBG calculation. Default value of 0.1.",
      false);

  std::stringstream ss;
  ss << "CalculateEWBG calculates the strength of the electroweak "
        "baryogenesis"
     << std::endl
     << "It is called either by " << std::endl
     << "./CalculateEWBG model input output FirstLine LastLine ConfigFile"
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
      arguments.emplace_back("--firstLine=" + std::string(argv[4]));
    }
    if (argc >= 6)
    {
      arguments.emplace_back("--lastLine=" + std::string(argv[5]));
    }
    if (argc >= 7)
    {
      arguments.emplace_back("--config=" + std::string(argv[6]));
    }
  }
  return arguments;
}
