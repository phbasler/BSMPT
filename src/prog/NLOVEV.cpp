// Copyright (C) 2020  Philipp Basler, Margarete Mühlleitner and Jonas Müller
// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

/**
 * @file
 * Calculates the VEV at T = 0 at NLO for a given Inputfile for a given subset
 * of lines in the file and adds it at the end of the line . One parameter point
 * per line.
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
using namespace std;
using namespace BSMPT;

struct CLIOptions
{
  BSMPT::ModelID::ModelIDs Model{ModelID::ModelIDs::NotSet};
  int FirstLine{}, LastLine{};
  std::string InputFile, OutputFile;
  bool UseGSL{Minimizer::UseGSLDefault};
  bool UseCMAES{Minimizer::UseLibCMAESDefault};
  bool UseNLopt{Minimizer::UseNLoptDefault};
  int WhichMinimizer{Minimizer::WhichMinimizerDefault};

  CLIOptions(int argc, char *argv[]);
  bool good() const;
};

int main(int argc, char *argv[])
try
{

  const CLIOptions args(argc, argv);
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

  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(args.Model);
  std::vector<double> Check;
  while (true)
  {

    getline(infile, linestr);
    if (linecounter > args.LastLine) break;
    if (linecounter == 1)
    {
      modelPointer->setUseIndexCol(linestr);
      outfile << linestr;
      auto legendCT = modelPointer->addLegendCT();
      for (auto x : legendCT)
        outfile << sep << x;
      auto legendVEV = modelPointer->addLegendVEV();
      for (auto x : legendVEV)
        outfile << sep << x;
      outfile << sep << "v_NLO";
      outfile << std::endl;
    }
    if (linecounter >= args.FirstLine and linecounter <= args.LastLine and
        linecounter != 1)
    {
      std::pair<std::vector<double>, std::vector<double>> parameters =
          modelPointer->initModel(linestr);

      if (args.FirstLine == args.LastLine) modelPointer->write();
      Check.clear();
      auto sol = Minimizer::Minimize_gen_all(modelPointer,
                                             0,
                                             Check,
                                             modelPointer->get_vevTreeMin(),
                                             args.WhichMinimizer);

      std::vector<double> solPot, solSym;
      solPot     = modelPointer->MinimizeOrderVEV(sol);
      double vev = modelPointer->EWSBVEV(solPot);

      outfile << linestr;
      outfile << sep << parameters.second << sep << sol << sep << vev
              << std::endl;

      if (args.FirstLine == args.LastLine)
      {
        auto dimensionnames = modelPointer->addLegendVEV();
        for (std::size_t i = 0; i < modelPointer->get_nVEV(); i++)
        {
          Logger::Write(LoggingLevel::Default,
                        dimensionnames.at(i) + " = " +
                            std::to_string(sol.at(i)) + " GeV");
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

CLIOptions::CLIOptions(int argc, char *argv[])
{
  std::vector<std::string> args;
  for (int i{1}; i < argc; ++i)
    args.push_back(argv[i]);

  if (argc < 6 or args.at(0) == "--help")
  {
    std::stringstream ss;
    int SizeOfFirstColumn = std::string("--TerminalOutput=           ").size();
    ss << "NLOVEV calculates the EW VEV at NLO" << std::endl
       << "It is called either by " << std::endl
       << "./NLOVEV model input output FirstLine LastLine" << std::endl
       << "or with the following arguments" << std::endl
       << std::setw(SizeOfFirstColumn) << std::left << "--help"
       << "Shows this menu" << std::endl
       << std::setw(SizeOfFirstColumn) << std::left << "--model="
       << "The model you want to investigate" << std::endl
       << std::setw(SizeOfFirstColumn) << std::left << "--input="
       << "The input file in tsv format" << std::endl
       << std::setw(SizeOfFirstColumn) << std::left << "--output="
       << "The output file in tsv format" << std::endl
       << std::setw(SizeOfFirstColumn) << std::left << "--FirstLine="
       << "The first line in the input file to calculate the NLO EW "
          "VEV. Expects line 1 to be a legend."
       << std::endl
       << std::setw(SizeOfFirstColumn) << std::left << "--LastLine="
       << "The last line in the input file to calculate the NLO EW VEV."
       << std::endl;
    std::string GSLhelp{"--UseGSL="};
    GSLhelp += Minimizer::UseGSLDefault ? "true" : "false";
    ss << std::setw(SizeOfFirstColumn) << std::left << GSLhelp
       << "Use the GSL library to minimize the effective potential"
       << std::endl;
    std::string CMAEShelp{"--UseCMAES="};
    CMAEShelp += Minimizer::UseLibCMAESDefault ? "true" : "false";
    ss << std::setw(SizeOfFirstColumn) << std::left << CMAEShelp
       << "Use the CMAES library to minimize the effective potential"
       << std::endl;
    std::string NLoptHelp{"--UseNLopt="};
    NLoptHelp += Minimizer::UseNLoptDefault ? "true" : "false";
    ss << std::setw(SizeOfFirstColumn) << std::left << NLoptHelp
       << "Use the NLopt library to minimize the effective potential"
       << std::endl;
    Logger::Write(LoggingLevel::Default, ss.str());
    ShowLoggerHelp();
    ShowInputError();
  }

  if (args.size() > 0 and args.at(0) == "--help")
  {
    throw int{0};
  }
  else if (argc < 6)
  {
    throw std::runtime_error("Too few arguments.");
  }

  const std::string prefix{"--"};
  bool UsePrefix = StringStartsWith(args.at(0), prefix);
  std::vector<std::string> UnusedArgs;
  if (UsePrefix)
  {
    for (const auto &arg : args)
    {
      auto el = arg;
      std::transform(el.begin(), el.end(), el.begin(), ::tolower);
      if (StringStartsWith(el, "--model="))
      {
        Model =
            BSMPT::ModelID::getModel(el.substr(std::string("--model=").size()));
      }
      else if (StringStartsWith(el, "--input="))
      {
        InputFile = arg.substr(std::string("--input=").size());
      }
      else if (StringStartsWith(el, "--output="))
      {
        OutputFile = arg.substr(std::string("--output=").size());
      }
      else if (StringStartsWith(el, "--firstline="))
      {
        FirstLine = std::stoi(el.substr(std::string("--firstline=").size()));
      }
      else if (StringStartsWith(el, "--lastline="))
      {
        LastLine = std::stoi(el.substr(std::string("--lastline=").size()));
      }
      else if (StringStartsWith(el, "--usegsl="))
      {
        UseGSL = el.substr(std::string("--usegsl=").size()) == "true";
      }
      else if (StringStartsWith(el, "--usecmaes="))
      {
        UseCMAES = el.substr(std::string("--usecmaes=").size()) == "true";
      }
      else if (StringStartsWith(el, "--usenlopt="))
      {
        UseNLopt = el.substr(std::string("--usenlopt=").size()) == "true";
      }
      else
      {
        UnusedArgs.push_back(el);
      }
    }
    WhichMinimizer = Minimizer::CalcWhichMinimizer(UseGSL, UseCMAES, UseNLopt);
    SetLogger(UnusedArgs);
  }
  else
  {
    Model      = ModelID::getModel(args.at(0));
    InputFile  = args.at(1);
    OutputFile = args.at(2);
    FirstLine  = std::stoi(args.at(3));
    LastLine   = std::stoi(args.at(4));
  }
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
