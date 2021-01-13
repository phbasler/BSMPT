/*
 * CreateMuGrid.cpp
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
 * Calculates the electroweak phase transition for a given Inputfile for a given
 * subset of lines in the file and adds it at the end of the line in the format
 * T_c v_c all single vevs. One parameter point per line.
 *
 */

#include "BSMPT/models/ClassPotentialOrigin.h" // for Class_Poten...
#include <BSMPT/baryo_calculation/transport_equations.h>
#include <BSMPT/minimizer/Minimizer.h>
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/utility.h>
#include <algorithm> // for copy, max
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>   // for shared_ptr
#include <stdlib.h> // for EXIT_FAILURE
#include <string>   // for getline
#include <utility>  // for pair
#include <vector>   // for vector
using namespace std;
using namespace BSMPT;

struct CLIOptions
{
  BSMPT::ModelID::ModelIDs Model{};
  int Line{};
  std::string InputFile, OutputFile;
  double vw{0.1};
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

  std::shared_ptr<Class_Potential_Origin> modelPointer =
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
  bool found      = false;
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
    return EXIT_FAILURE;
  }

  std::vector<double> parStart, parEnd;
  parStart = std::vector<double>(modelPointer->get_NHiggs(), 0);
  auto EWPT =
      Minimizer::PTFinder_gen_all(modelPointer, 0, 300, args.WhichMinimizer);

  // find the minimum in the symmetric phase. For this minimise at T = Tc + 1
  std::vector<double> vevsymmetricSolution, checksym, startpoint;
  for (const auto &el : EWPT.EWMinimum)
    startpoint.push_back(0.5 * el);
  vevsymmetricSolution = Minimizer::Minimize_gen_all(
      modelPointer, EWPT.Tc + 1, checksym, startpoint, args.WhichMinimizer);

  double absvevsymmetricSolution = 0;
  for (const auto &x : vevsymmetricSolution)
    absvevsymmetricSolution += std::pow(x, 2);

  if (absvevsymmetricSolution != 0)
  {
    std::cout
        << "Check for sign of non vanishing components of symmetric vacuum : "
        << std::endl;
  }

  struct Baryo::GSL_integration_mubl p;
  p.init(args.vw,
         EWPT.EWMinimum,
         vevsymmetricSolution,
         EWPT.Tc,
         modelPointer,
         args.WhichMinimizer);

  std::cout << "vw = " << args.vw << std::endl;
  std::cout << "LW = " << p.getLW() * EWPT.Tc << "/TC" << std::endl;
  std::cout << "T_C = " << EWPT.Tc << std::endl;
  for (std::size_t i = 0; i < modelPointer->get_nVEV(); i++)
    std::cout << "v_" << i << " = " << EWPT.EWMinimum.at(i) << std::endl;

  std::size_t nstep = 1000;
  double zmin       = 0;
  double stepsize   = (p.getZMAX() - zmin) / nstep;
  outfile << "z\tmu_{B_L}" << std::endl;
  for (std::size_t i = 0; i <= nstep; i++)
  {
    double z = zmin + stepsize * i;
    outfile << z << sep << Baryo::mubl_func(z, &p) << std::endl;
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

  if (argc < 5 or args.at(0) == "--help")
  {
    int SizeOfFirstColumn = std::string("--TerminalOutput=           ").size();
    std::cout << "CreateMuGrid calculates the mu_{BL} potential in front of "
                 "the bubble wall"
              << std::endl
              << "It is called either by " << std::endl
              << "./CreateMuGrid model input output Line" << std::endl
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
              << "The line in the input file to calculate the EWPT. Expects "
                 "line 1 to be a legend."
              << std::endl
              << std::setw(SizeOfFirstColumn) << std::left << "--vw="
              << "Wall velocity for the EWBG calculation. Default value of 0.1."
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
  else if (argc < 5)
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
      else if (StringStartsWith(el, "--vw="))
      {
        vw = std::stod(el.substr(std::string("--vw=").size()));
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
    Model      = ModelID::getModel(args.at(0));
    InputFile  = args.at(1);
    OutputFile = args.at(2);
    Line       = std::stoi(args.at(3));
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

  if (vw <= 0 or vw > 1)
  {
    throw std::runtime_error("The wall velocity has to be between 0 and 1.");
  }
  if (Model == ModelID::ModelIDs::NotSet)
  {

    std::cerr
        << "Your Model parameter does not match with the implemented Models."
        << std::endl;
    ShowInputError();
    return false;
  }
  if (Line < 1)
  {
    std::cerr << "Start line counting with 1" << std::endl;
    return false;
  }
  return true;
}
