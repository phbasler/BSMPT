/*
 * VEVEVO.cpp
 *
 *
 *      Copyright (C) 2020  Philipp Basler, Margarete Mühlleitner and Jonas
 Müller

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * @file
 * This program calculates the development of the VeVs with the Temperature.
 */

#include <BSMPT/minimizer/Minimizer.h>
#include <BSMPT/models/ClassPotentialOrigin.h> // for Class_Potential_Origin
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/utility.h>
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

  std::vector<double> sol, start, solPot;

  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(args.Model);
  std::ifstream infile(args.InputFile);
  if (!infile.good())
  {
    std::cout << "Input file not found " << std::endl;
    return EXIT_FAILURE;
  }
  std::ofstream outfile(args.OutputFile);
  if (!outfile.good())
  {
    std::cout << "Can not create file " << args.OutputFile << std::endl;
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
    std::cout << "Line not found !\n";
    return -1;
  }

  std::vector<double> Check;
  double vev{0.0};

  std::cout << std::scientific;
  std::cout << std::setprecision(16);
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
    sol = Minimizer::Minimize_gen_all(
        modelPointer, Temp, Check, start, args.WhichMinimizer);
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
  std::cerr << e.what() << std::endl;
  return EXIT_FAILURE;
}

CLIOptions::CLIOptions(int argc, char *argv[])
{
  std::vector<std::string> args;
  for (int i{1}; i < argc; ++i)
    args.push_back(argv[i]);

  if (argc < 8 or args.at(0) == "--help")
  {
    int SizeOfFirstColumn =
        std::string("--TemperatureStart=           ").size();
    std::cout << "VEVEVO calculates the evolution of the global minimum with "
                 "rising temperature for a given parameter point"
              << std::endl
              << "It is called either by " << std::endl
              << "./VEVEVO Model Inputfile Outputfile Line TemperatureStart "
                 "TemperatureStep TemperatureEnd"
              << std::endl
              << "or with the following arguments" << std::endl
              << std::setw(SizeOfFirstColumn) << std::left << "--help"
              << "Shows this menu" << std::endl
              << std::setw(SizeOfFirstColumn) << std::left << "--model="
              << "The model you want to investigate" << std::endl
              << std::setw(SizeOfFirstColumn) << std::left << "--input="
              << "The input file in tsv format" << std::endl
              << std::setw(SizeOfFirstColumn) << std::left << "--output="
              << "The output file in tsv format" << std::endl
              << std::setw(SizeOfFirstColumn) << std::left << "--Line="
              << "The line in the input file with the given parameter point. "
                 "Expects line 1 to be a legend."
              << std::endl
              << std::setw(SizeOfFirstColumn) << std::left
              << "--TemperatureStart="
              << "The starting temperature to calculate the global minimum."
              << std::endl
              << std::setw(SizeOfFirstColumn) << std::left
              << "--TemperatureStep="
              << "The stepsize for the temperature." << std::endl
              << std::setw(SizeOfFirstColumn) << std::left
              << "--TemperatureEnd="
              << "The last temperature to calculate the global minimum."
              << std::endl;
    std::string GSLhelp{"--UseGSL="};
    GSLhelp += Minimizer::UseGSLDefault ? "true" : "false";
    std::cout << std::setw(SizeOfFirstColumn) << std::left << GSLhelp
              << "Use the GSL library to minimize the effective potential"
              << std::endl;
    std::string CMAEShelp{"--UseCMAES="};
    CMAEShelp += Minimizer::UseLibCMAESDefault ? "true" : "false";
    std::cout << std::setw(SizeOfFirstColumn) << std::left << CMAEShelp
              << "Use the CMAES library to minimize the effective potential"
              << std::endl;
    std::string NLoptHelp{"--UseNLopt="};
    NLoptHelp += Minimizer::UseNLoptDefault ? "true" : "false";
    std::cout << std::setw(SizeOfFirstColumn) << std::left << NLoptHelp
              << "Use the NLopt library to minimize the effective potential"
              << std::endl;
    ShowInputError();
  }

  if (args.size() > 0 and args.at(0) == "--help")
  {
    throw int{0};
  }
  else if (argc < 8)
  {
    throw std::runtime_error("Too few arguments.");
  }

  const std::string prefix{"--"};
  bool UsePrefix = StringStartsWith(args.at(0), prefix);
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
      else if (StringStartsWith(el, "--line="))
      {
        Line = std::stoi(el.substr(std::string("--line=").size()));
      }
      else if (StringStartsWith(el, "--temperaturestart="))
      {
        TemperatureStart =
            std::stod(el.substr(std::string("--temperaturestart=").size()));
      }
      else if (StringStartsWith(el, "--temperaturestep="))
      {
        TemperatureStep =
            std::stod(el.substr(std::string("--temperaturestep=").size()));
      }
      else if (StringStartsWith(el, "--temperatureend="))
      {
        TemperatureEnd =
            std::stod(el.substr(std::string("--temperatureend=").size()));
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
    }
    WhichMinimizer = Minimizer::CalcWhichMinimizer(UseGSL, UseCMAES, UseNLopt);
  }
  else
  {
    Model            = ModelID::getModel(args.at(0));
    InputFile        = args.at(1);
    OutputFile       = args.at(2);
    Line             = std::stoi(args.at(3));
    TemperatureStart = std::stod(args.at(4));
    TemperatureStep  = std::stod(args.at(5));
    TemperatureEnd   = std::stod(args.at(6));
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

    std::cerr
        << "Your Model parameter does not match with the implemented Models."
        << std::endl;
    ShowInputError();
    return false;
  }
  return true;
}
