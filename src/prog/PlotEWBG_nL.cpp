// Copyright (C) 2020  Philipp Basler, Margarete Mühlleitner and Jonas
// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

/**
 * @file
 * This program calculates the left-handed density in front of the bubble wall
 * as a function of the distance z in both approaches.
 */

#include <BSMPT/baryo_calculation/CalculateEtaInterface.h>
#include <BSMPT/baryo_calculation/Fluid_Type/tau_source.h>
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
using namespace Baryo;

struct CLIOptions
{
  BSMPT::ModelID::ModelIDs Model{ModelID::ModelIDs::NotSet};
  int Line{};
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
      ModelID::FChoose(args.Model);

  std::vector<double> start, solPot;

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
      outfile << "T_c_nLvar" << sep << "omega_c_nLvar" << sep << "vw_nLvar"
              << sep << "LW_nLvar" << sep;
      outfile << "z" << sep << "z/LW" << sep << "nL_VIA" << sep << "muL_FH"
              << sep << "nL_FH";
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
    Logger::Write(
        LoggingLevel::ProgDetailed, "PTFinder called...", __FILE__, __LINE__);
  auto EWPT =
      Minimizer::PTFinder_gen_all(modelPointer, 0, 300, args.WhichMinimizer);
  // SFOEWPT FOUND
  if (EWPT.StatusFlag == Minimizer::MinimizerStatus::SUCCESS and
      C_PT * EWPT.Tc < EWPT.vc)
  {
    if (args.TerminalOutput)
      Logger::Write(
          LoggingLevel::ProgDetailed, "SFOEWPT found...", __FILE__, __LINE__);
    std::vector<double> vcritical, vbarrier;
    vcritical = EWPT.EWMinimum;
    double TC = EWPT.Tc;
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

    /////////////////////////////////////////////////////////////////////////////////
    std::size_t nstep = 100;

    if (args.TerminalOutput)
      Logger::Write(LoggingLevel::ProgDetailed,
                    "Set up the numerics ",
                    __FILE__,
                    __LINE__);
    EtaInterface.setNumerics(
        args.vw,
        EWPT.EWMinimum,
        vevsymmetricSolution,
        TC,
        modelPointer,
        args.WhichMinimizer); // Set up parameter container for Baryo
                              // Calculation-->Calls container.init
    if (args.TerminalOutput)
      Logger::Write(LoggingLevel::ProgDetailed,
                    "Starting setting the class instances",
                    __FILE__,
                    __LINE__);
    BSMPT::Baryo::tau_source C_tau(SMConstants);
    EtaInterface.set_transport_method(
        TransportMethod::tau); // setting to tau class
    bool botflag     = true;
    auto class_GamM  = EtaInterface.get_class_CalcGamM();
    auto class_ScP   = EtaInterface.get_class_Scp();
    auto class_kappa = EtaInterface.get_class_kappa();
    auto Integration_mubl{EtaInterface.getGSL_integration_mubl_container()};
    C_tau.set_class(
        botflag, Integration_mubl, class_GamM, class_ScP, class_kappa);

    auto tau_arr_nL = set_up_nL_grid(nstep, Integration_mubl, C_tau);
    auto FH_GSL =
        generate_mubl_spline(Integration_mubl, static_cast<int>(nstep));

    ///////
    /// Outfile
    ///////
    for (std::size_t i = 0; i < nstep; i++)
    {
      outfile << linestr << sep;
      outfile << TC << sep << EWPT.vc << sep << args.vw << sep
              << EtaInterface.getLW() << sep;
      outfile << tau_arr_nL.first.at(i) << sep;
      outfile << tau_arr_nL.first.at(i) / EtaInterface.getLW() << sep;
      outfile << tau_arr_nL.second.at(i) << sep;
      outfile << FH_GSL.spline(tau_arr_nL.first.at(i)) << sep;
      outfile << FH_GSL.spline(tau_arr_nL.first.at(i)) * std::pow(TC, 2);
      outfile << std::endl;
    }
    /////////////////////////////////////////////////////////////////////////////////

  } // END: SFOEWPT FOUND
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
  if (Line < 1)
  {
    Logger::Write(LoggingLevel::Default, "Start line counting with 1");
    return false;
  }
  return true;
}

CLIOptions::CLIOptions(const BSMPT::parser &argparser)
{
  argparser.check_required_parameters();
  Model      = BSMPT::ModelID::getModel(argparser.get_value("model"));
  InputFile  = argparser.get_value("input");
  OutputFile = argparser.get_value("output");
  Line       = argparser.get_value<int>("line");
  vw         = argparser.get_value<double>("vw");
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
      "vw",
      "Wall velocity for the EWBG calculation. Default value of 0.1.",
      true);
  argparser.add_argument("config", "The EWBG config file.", true);
  argparser.add_argument(
      "terminalOutput",
      "y/n Turns on additional information in the terminal during "
      "the calculation.",
      false);

  std::stringstream ss;
  ss << "Calculation of the left-handed chemical potentials or "
        "particle densities triggering the EW sphaleron transitions as a "
        "function of the wall distance z ."
     << std::endl
     << "It is called either by " << std::endl
     << "./PlotEWBG_nL Model Inputfile Outputfile Line vw "
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
      arguments.emplace_back("--vw=" + std::string(argv[5]));
    }
    if (argc >= 7)
    {
      arguments.emplace_back("--config=" + std::string(argv[6]));
    }
  }
  return arguments;
}
