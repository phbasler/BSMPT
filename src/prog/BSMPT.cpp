// Copyright (C) 2020  Philipp Basler, Margarete Mühlleitner and Jonas
// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

/**
 * @file
 * Calculates the electroweak phase transition for a given Inputfile for a given
 * subset of lines in the file and adds it at the end of the line in the format
 * T_c v_c all single vevs. One parameter point per line.
 *
 */

#include <BSMPT/minimizer/Minimizer.h>
#include <BSMPT/models/ClassPotentialOrigin.h> // for Class_Potential_Origin
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/utility/Logger.h>
#include <BSMPT/utility/utility.h>
#include <algorithm> // for copy, max
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>   // for shared_ptr, __shared_...
#include <stdlib.h> // for atoi, EXIT_FAILURE
#include <string>   // for string, operator<<
#include <utility>  // for pair
#include <vector>   // for vector

#include <BSMPT/utility/parser.h>

using namespace std;
using namespace BSMPT;

struct CLIOptions
{
  BSMPT::ModelID::ModelIDs Model{ModelID::ModelIDs::NotSet};
  int FirstLine{0}, LastLine{0};
  std::string InputFile, OutputFile;
  bool TerminalOutput{false};
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
  /**
   * PrintErrorLines decides if parameter points with no valid EWPT (no NLO
   * stability or T=300 vanishing VEV) are printed in the output file
   */
  bool PrintErrorLines = true;

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
    Logger::Write(LoggingLevel::Default,
                  "Input file " + args.InputFile + " not found ");
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
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(args.Model);
  while (getline(infile, linestr))
  {
    if (linecounter > args.LastLine) break;

    if (linecounter == 1)
    {
      outfile << linestr << sep << modelPointer->addLegendCT() << sep
              << modelPointer->addLegendTemp() << std::endl;

      modelPointer->setUseIndexCol(linestr);
    }
    if (linecounter >= args.FirstLine and linecounter <= args.LastLine and
        linecounter != 1)
    {
      if (args.TerminalOutput)
      {
        Logger::Write(LoggingLevel::ProgDetailed,
                      "Currently at line " + std::to_string(linecounter));
      }
      std::pair<std::vector<double>, std::vector<double>> parameters =
          modelPointer->initModel(linestr);
      if (args.FirstLine == args.LastLine)
      {
        modelPointer->write();
      }

      auto EWPT = Minimizer::PTFinder_gen_all(
          modelPointer, 0, 300, args.WhichMinimizer, args.UseMultithreading);
      std::vector<double> vevsymmetricSolution, checksym, startpoint;
      for (const auto &el : EWPT.EWMinimum)
        startpoint.push_back(0.5 * el);
      auto VEVsym = Minimizer::Minimize_gen_all(modelPointer,
                                                EWPT.Tc + 1,
                                                checksym,
                                                startpoint,
                                                args.WhichMinimizer,
                                                args.UseMultithreading);

      if (args.FirstLine == args.LastLine)
      {
        auto dimensionnames = modelPointer->addLegendTemp();
        Logger::Write(
            LoggingLevel::Default,
            "Success ? " + std::to_string(static_cast<int>(EWPT.StatusFlag)) +
                sep + " (1 = Yes , -1 = No, v/T reached a value below " +
                std::to_string(C_PT) + " during the calculation)");
        if (EWPT.StatusFlag == Minimizer::MinimizerStatus::SUCCESS)
        {
          Logger::Write(LoggingLevel::Default,
                        dimensionnames.at(1) + " = " + std::to_string(EWPT.vc) +
                            " GeV");
          Logger::Write(LoggingLevel::Default,
                        dimensionnames.at(0) + " = " + std::to_string(EWPT.Tc) +
                            " GeV");
          Logger::Write(LoggingLevel::Default,
                        "xi_c = " + dimensionnames.at(2) + " = " +
                            std::to_string(EWPT.vc / EWPT.Tc));
          for (std::size_t i = 0; i < modelPointer->get_nVEV(); i++)
          {
            Logger::Write(LoggingLevel::Default,
                          dimensionnames.at(i + 3) + " = " +
                              std::to_string(EWPT.EWMinimum.at(i)) + " GeV");
          }
          Logger::Write(LoggingLevel::Default, "Symmetric VEV config");
          for (std::size_t i = 0; i < modelPointer->get_nVEV(); i++)
          {
            Logger::Write(LoggingLevel::Default,
                          dimensionnames.at(i + 3) + " = " +
                              std::to_string(VEVsym.at(i)) + " GeV");
          }
        }
        else if (EWPT.StatusFlag ==
                 Minimizer::MinimizerStatus::NOTVANISHINGATFINALTEMP)
        {
          Logger::Write(
              LoggingLevel::Default,
              dimensionnames.at(1) +
                  " != 0 GeV at T = 300 GeV. No SFOEWPT is possible.");
        }
        else if (EWPT.StatusFlag == Minimizer::MinimizerStatus::NLOVEVZEROORINF)
        {
          Logger::Write(LoggingLevel::Default,
                        dimensionnames.at(1) + " = 0 / > 255 GeV at T = 0 GeV. "
                                               "The point is not NLO stable.");
        }
        else if (EWPT.StatusFlag == Minimizer::MinimizerStatus::NOTNLOSTABLE)
        {
          Logger::Write(
              LoggingLevel::Default,
              dimensionnames.at(1) +
                  " != vEW GeV at T = 0 GeV. The point is not NLO stable.");
        }
        else if (EWPT.StatusFlag ==
                 Minimizer::MinimizerStatus::NUMERICALLYUNSTABLE)
        {
          Logger::Write(LoggingLevel::Default,
                        "The point is numerically unstable.");
        }
        else if (EWPT.StatusFlag == Minimizer::MinimizerStatus::BELOWTHRESHOLD)
        {
          Logger::Write(LoggingLevel::Default,
                        dimensionnames.at(1) + " < " + std::to_string(C_PT) +
                            " found.");
        }
      }
      if (PrintErrorLines)
      {
        outfile << linestr;
        outfile << sep << parameters.second;
        outfile << sep << EWPT.Tc << sep << EWPT.vc;
        if (EWPT.vc > C_PT * EWPT.Tc and
            EWPT.StatusFlag == Minimizer::MinimizerStatus::SUCCESS)
          outfile << sep << EWPT.vc / EWPT.Tc;
        else
          outfile << sep << static_cast<int>(EWPT.StatusFlag);
        outfile << sep << EWPT.EWMinimum;
        outfile << std::endl;
      }
      else if (EWPT.StatusFlag == Minimizer::MinimizerStatus::SUCCESS)
      {
        if (C_PT * EWPT.Tc < EWPT.vc)
        {
          outfile << linestr << sep << parameters.second;
          outfile << sep << EWPT.Tc << sep << EWPT.vc;
          outfile << sep << EWPT.vc / EWPT.Tc;
          outfile << sep << EWPT.EWMinimum;
          outfile << std::endl;
        }
      }
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

  try
  {
    UseGSL = argparser.get_value<bool>("UseGSL");
  }
  catch (BSMPT::parserException &)
  {
  }

  try
  {
    UseCMAES = argparser.get_value<bool>("UseCMAES");
  }
  catch (BSMPT::parserException &)
  {
  }

  try
  {
    UseNLopt = argparser.get_value<bool>("UseNLopt");
  }
  catch (BSMPT::parserException &)
  {
  }

  try
  {
    UseMultithreading = argparser.get_value<bool>("UseMultithreading");
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
  ss << "BSMPT calculates the strength of the electroweak phase transition"
     << std::endl
     << "It is called either by " << std::endl
     << "./BSMPT model input output FirstLine LastLine" << std::endl
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
      arguments.emplace_back("--FirstLine=" + std::string(argv[4]));
    }
    if (argc >= 6)
    {
      arguments.emplace_back("--LastLine=" + std::string(argv[5]));
    }
  }
  return arguments;
}
