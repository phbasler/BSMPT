// Copyright (C) 2024 Lisa Biermann, Margarete Mühlleitner, Rui Santos, João
// Viana SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and
// Jonas Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

/**
 * @file
 * This program traces all minima in a temperature range
 *
 */

#include <BSMPT/minimizer/Minimizer.h>
#include <BSMPT/minimum_tracer/minimum_tracer.h> // MinimumTracer
#include <BSMPT/models/ClassPotentialOrigin.h>   // for Class_Potential_Origin
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/transition_tracer/transition_tracer.h> // TransitionTracer
#include <BSMPT/utility/Logger.h>
#include <BSMPT/utility/parser.h>
#include <BSMPT/utility/utility.h>
#include <algorithm> // for copy, max
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory> // for shared_ptr, __shared_...
#include <random>
#include <stdlib.h> // for atoi, EXIT_FAILURE
#include <string>   // for string, operator<<
#include <utility>  // for pair
#include <vector>   // for vector

using namespace std;
using namespace BSMPT;

struct CLIOptions
{
  BSMPT::ModelID::ModelIDs Model{ModelID::ModelIDs::NotSet};
  int line{0}, npoints{100};
  std::string inputfile, outputfile;
  double templow{0}, temphigh{300};
  bool UseGSL{Minimizer::UseGSLDefault};
  bool UseCMAES{Minimizer::UseLibCMAESDefault};
  bool UseNLopt{Minimizer::UseNLoptDefault};
  int WhichMinimizer{Minimizer::WhichMinimizerDefault};
  bool UseMultithreading{false};
  int UseMultiStepPTMode{-1};
  int CheckEWSymmetryRestoration{1};
  int num_check_pts{10};

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

  std::ifstream infile(args.inputfile);
  if (!infile.good())
  {
    Logger::Write(LoggingLevel::Default,
                  "Input file " + args.inputfile + " not found ");
    return EXIT_FAILURE;
  }

  Logger::Write(LoggingLevel::ProgDetailed, "Found file");

  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(args.Model, SMConstants);

  Logger::Write(LoggingLevel::ProgDetailed, "Created modelpointer ");

  std::string linestr, linestr_store;
  int linecounter = 1, filecounter = 1;

  while (getline(infile, linestr))
  {
    if (linecounter == 1) linestr_store = linestr;

    if (linecounter > args.line) break;

    if (linecounter == args.line)
    {
      Logger::Write(LoggingLevel::ProgDetailed,
                    "Currently at line " + std::to_string(linecounter));

      std::string outfilename =
          args.outputfile + "_" + std::to_string(filecounter) + ".tsv";
      Logger::Write(LoggingLevel::ProgDetailed,
                    "Creating outfile with name " + outfilename);
      std::ofstream outfile(outfilename);
      if (!outfile.good())
      {
        Logger::Write(LoggingLevel::Default,
                      "Can not create file " + outfilename);
        return EXIT_FAILURE;
      }

      modelPointer->setUseIndexCol(linestr_store);

      std::pair<std::vector<double>, std::vector<double>> parameters =
          modelPointer->initModel(linestr);

      modelPointer->write();

      auto start = std::chrono::high_resolution_clock::now();

      std::shared_ptr<MinimumTracer> MinTracer(new MinimumTracer(
          modelPointer, args.WhichMinimizer, args.UseMultithreading));

      // NLO stability check
      bool nlostable = modelPointer->CheckNLOVEV(
          MinTracer->ConvertToVEVDim(MinTracer->GetGlobalMinimum(0)));
      StatusNLOStability status_nlostable =
          MinTracer->GetStatusNLOVEV(nlostable);
      Logger::Write(LoggingLevel::ProgDetailed,
                    "Status of NLO stability check is: " +
                        StatusNLOStabilityToString.at(status_nlostable));

      // EWSR check
      double EWSymmetryRestoration_status = 0;
      StatusEWSR status_ewsr              = StatusEWSR::Off;

      if (args.CheckEWSymmetryRestoration > 0)
      {
        EWSymmetryRestoration_status =
            MinTracer->IsThereEWSymmetryRestoration();
        status_ewsr = MinTracer->GetStatusEWSR(EWSymmetryRestoration_status);
      }
      else
      {
        Logger::Write(LoggingLevel::ProgDetailed,
                      "Check for EW symmetry restoration is disabled.\n");
      }

      // phase tracking
      Logger::Write(
          LoggingLevel::ProgDetailed,
          "Track phases in between T_low = " + std::to_string(args.templow) +
              " and T_high = " + std::to_string(args.temphigh));

      Vacuum vac(args.templow,
                 args.temphigh,
                 MinTracer,
                 modelPointer,
                 args.UseMultiStepPTMode,
                 args.num_check_pts);

      Logger::Write(LoggingLevel::ProgDetailed,
                    "Found and traced " +
                        std::to_string(vac.PhasesList.size()) +
                        " minima with status = " +
                        StatusTracingToString.at(vac.status_vacuum) +
                        ".\n-------------------------------");

      if (vac.PhasesList.size() != 2)
        throw std::runtime_error(
            "This code is not ready for handling anything other than "
            "two phases!\n");

      // prepare legend
      std::vector<std::string> LegendMinima;
      LegendMinima.push_back("status_nlo_stability");
      LegendMinima.push_back("status_ewsr");
      LegendMinima.push_back("status_tracing");

      LegendMinima.push_back("T");
      for (std::size_t j = 0; j < modelPointer->get_nVEV(); j++)
        LegendMinima.push_back(modelPointer->addLegendVEV().at(j));
      for (std::size_t j = 0; j < modelPointer->get_NHiggs(); j++)
        LegendMinima.push_back("mS2_" + to_string(j));
      for (std::size_t j = 0; j < modelPointer->get_NLepton(); j++)
        LegendMinima.push_back("mL2_" + to_string(j));
      for (std::size_t j = 0; j < modelPointer->get_NQuarks(); j++)
        LegendMinima.push_back("mQ2_" + to_string(j));
      for (std::size_t j = 0; j < modelPointer->get_NGauge(); j++)
        LegendMinima.push_back("mG2_" + to_string(j));

      LegendMinima.push_back("runtime");
      outfile << linestr_store << sep << modelPointer->addLegendCT() << sep
              << LegendMinima << std::endl;

      std::size_t length = 0;
      for (std::size_t j = 0; j < vac.PhasesList.size(); j++)
      {
        if (length == 0)
        {
          length = vac.PhasesList.at(j).MinimumPhaseVector.size();
        }
        else if (length < vac.PhasesList.at(j).MinimumPhaseVector.size())
        {
          length = vac.PhasesList.at(j).MinimumPhaseVector.size();
        }
      }

      if (length == 0)
      {
        outfile << std::setprecision(16);
        outfile << linestr;
        outfile << sep << parameters.second;
        outfile << sep << status_nlostable;
        outfile << sep << status_ewsr;
        outfile << sep << vac.status_vacuum;
        auto time = std::chrono::duration_cast<std::chrono::milliseconds>(
                        std::chrono::high_resolution_clock::now() - start)
                        .count() /
                    1000.;

        outfile << sep << time;
        outfile << std::endl;
      }

      std::vector<double> T_list;
      for (int n = 0; n < args.npoints; n++)
        T_list.push_back(args.templow +
                         n * (args.temphigh - args.templow) / args.npoints);

      T_list.insert(std::lower_bound(T_list.begin(),
                                     T_list.end(),
                                     vac.CoexPhasesList.at(0).crit_temp),
                    vac.CoexPhasesList.at(0).crit_temp);

      for (const auto T : T_list)
      {
        outfile << std::setprecision(16);
        outfile << linestr;
        outfile << sep << parameters.second;
        outfile << sep << status_nlostable;
        outfile << sep << status_ewsr;
        outfile << sep << vac.status_vacuum;
        outfile << sep << T;

        std::vector<double> vev;

        if (T <= vac.CoexPhasesList.at(0).crit_temp)
        {
          vev = vac.PhasesList.at(1).Get(T).point;
        }
        else
        {
          vev = vac.PhasesList.at(0).Get(T).point;
        }

        outfile << sep << vev;
        outfile << sep
                << modelPointer->HiggsMassesSquared(
                       modelPointer->MinimizeOrderVEV(vev), T);
        outfile << sep
                << modelPointer->LeptonMassesSquared(
                       modelPointer->MinimizeOrderVEV(vev));
        outfile << sep
                << modelPointer->QuarkMassesSquared(
                       modelPointer->MinimizeOrderVEV(vev));
        outfile << sep
                << modelPointer->GaugeMassesSquared(
                       modelPointer->MinimizeOrderVEV(vev), T);

        auto time = std::chrono::duration_cast<std::chrono::milliseconds>(
                        std::chrono::high_resolution_clock::now() - start)
                        .count() /
                    1000.;

        outfile << sep << time;
        outfile << std::endl;
      }

      filecounter++;
      outfile.close();

      auto time = std::chrono::duration_cast<std::chrono::milliseconds>(
                      std::chrono::high_resolution_clock::now() - start)
                      .count() /
                  1000.;

      BSMPT::Logger::Write(BSMPT::LoggingLevel::ProgDetailed,
                           "\nTook\t" + std::to_string(time) + " seconds.\n");
    }
    linecounter++;
    if (infile.eof()) break;
  }
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

  if (Model == ModelID::ModelIDs::NotSet)
  {
    Logger::Write(
        LoggingLevel::Default,
        "Your Model parameter does not match with the implemented Models.");
    ShowInputError();
    return false;
  }
  if (line == 0)
  {
    Logger::Write(LoggingLevel::Default, "line not set.");
    return false;
  }
  if (templow >= temphigh)
  {
    Logger::Write(LoggingLevel::Default,
                  "Invalid temperature choice. Thigh has to be > 0 GeV.");
    return false;
  }
  if (npoints <= 0)
  {
    Logger::Write(LoggingLevel::Default,
                  "Invalid npoints choice. npoints has to be > 0.");
    return false;
  }
  if (UseMultiStepPTMode > 3 or UseMultiStepPTMode < -1)
  {
    Logger::Write(LoggingLevel::Default, "Invalid choice for MultiStepPTMode.");
    return false;
  }
  if (num_check_pts < 0)
  {
    Logger::Write(LoggingLevel::Default, "Invalid choice for num_check_pts.");
    return false;
  }
  if (CheckEWSymmetryRestoration > 2 or CheckEWSymmetryRestoration < 0)
  {
    Logger::Write(LoggingLevel::Default,
                  "Invalid choice for CheckEWSymmetryRestoration.");
    return false;
  }

  return true;
}

CLIOptions::CLIOptions(const BSMPT::parser &argparser)
{
  std::stringstream ss;
  argparser.check_required_parameters();

  // required arguments
  Model      = BSMPT::ModelID::getModel(argparser.get_value("model"));
  inputfile  = argparser.get_value("input");
  outputfile = argparser.get_value("output");
  line       = argparser.get_value<int>("line");

  // optional arguments
  try
  {
    temphigh = argparser.get_value<double>("thigh");
  }
  catch (BSMPT::parserException &)
  {
    ss << "--thigh not set, using default value: " << temphigh << "\n";
  }

  try
  {
    npoints = argparser.get_value<double>("npoints");
  }
  catch (BSMPT::parserException &)
  {
    ss << "--npoints not set, using default value: " << temphigh << "\n";
  }

  std::string GSLhelp   = Minimizer::UseGSLDefault ? "true" : "false";
  std::string CMAEShelp = Minimizer::UseLibCMAESDefault ? "true" : "false";
  std::string NLoptHelp = Minimizer::UseNLoptDefault ? "true" : "false";
  try
  {
    UseGSL = (argparser.get_value("usegsl") == "true");
  }
  catch (BSMPT::parserException &)
  {
    ss << "--usegsl not set, using default value: " << GSLhelp << "\n";
  }

  try
  {
    UseCMAES = (argparser.get_value("usecmaes") == "true");
  }
  catch (BSMPT::parserException &)
  {
    ss << "--usecmaes not set, using default value: " << CMAEShelp << "\n";
  }

  try
  {
    UseNLopt = (argparser.get_value("usenlopt") == "true");
  }
  catch (BSMPT::parserException &)
  {
    ss << "--usenlopt not set, using default value: " << NLoptHelp << "\n";
  }

  try
  {
    UseMultithreading = (argparser.get_value("usemultithreading") == "true");
  }
  catch (BSMPT::parserException &)
  {
    ss << "--usemultithreading not set, using default value: false\n";
  }

  // UseMultiStepPTMode
  try
  {
    auto multistepPT_string = argparser.get_value("multistepmode");
    if (multistepPT_string == "default")
    {
      UseMultiStepPTMode = -1;
    }
    else if (multistepPT_string == "auto")
    {
      UseMultiStepPTMode = 3;
    }
    else
    {
      UseMultiStepPTMode = std::stoi(multistepPT_string);
    }
  }
  catch (BSMPT::parserException &)
  {
    ss << "--multistepmode not set, using default value: default\n";
  }

  try
  {
    num_check_pts = argparser.get_value<int>("num_pts");
  }
  catch (BSMPT::parserException &)
  {
    ss << "--num_check_pts not set, using default value: " << num_check_pts
       << "\n";
  }

  // CheckEWSymmetryRestoration
  try
  {
    auto ewsr_string = argparser.get_value("checkewsr");
    if (ewsr_string == "on")
    {
      CheckEWSymmetryRestoration = 1;
    }
    else if (ewsr_string == "off")
    {
      CheckEWSymmetryRestoration = 0;
    }
  }
  catch (BSMPT::parserException &)
  {
    ss << "--checkewsr not set, using default value: on\n";
  }

  WhichMinimizer = Minimizer::CalcWhichMinimizer(UseGSL, UseCMAES, UseNLopt);

  Logger::Write(LoggingLevel::ProgDetailed, ss.str());
}

BSMPT::parser prepare_parser()
{
  BSMPT::parser argparser(true);
  argparser.add_argument("help", "shows this menu", false);
  argparser.add_argument("model", "[*] model name", true);
  argparser.add_argument("input", "[*] input file (in tsv format)", true);
  argparser.add_argument("output", "[*] output file (in tsv format)", true);
  argparser.add_argument(
      "line", "[*] line number the parameter in input file", true);
  argparser.add_subtext("(expects line 1 to be a legend)");
  argparser.add_argument("npoints", "number temperature points", "100", false);
  argparser.add_argument("thigh", "high temperature [GeV]", "300", false);
  argparser.add_argument(
      "multistepmode", "multi-step PT mode", "default", false);
  argparser.add_subtext("default: default mode");
  argparser.add_subtext("0: single-step PT mode");
  argparser.add_subtext(">0 for multi-step PT modes:");
  argparser.add_subtext("1: tracing coverage");
  argparser.add_subtext("2: global minimum tracing coverage");
  argparser.add_subtext("auto: automatic mode");
  argparser.add_argument(
      "num_pts", "intermediate grid-size for default mode", "10", false);
  argparser.add_argument(
      "checkewsr", "check for EWSR at high temperature", "on", false);
  argparser.add_subtext("on: perform check");
  argparser.add_subtext("off: check disabled");

  std::string GSLhelp   = Minimizer::UseGSLDefault ? "true" : "false";
  std::string CMAEShelp = Minimizer::UseLibCMAESDefault ? "true" : "false";
  std::string NLoptHelp = Minimizer::UseNLoptDefault ? "true" : "false";

  argparser.add_argument(
      "usegsl", "use GSL library for minimization", GSLhelp, false);
  argparser.add_argument(
      "usecmaes", "use CMAES library  for minimization", CMAEShelp, false);
  argparser.add_argument(
      "usenlopt", "use NLopt library for minimization", NLoptHelp, false);
  argparser.add_argument("usemultithreading",
                         "enable multi-threading for minimizers",
                         "false",
                         false);
  argparser.add_argument(
      "json", "use a json file instead of cli parameters", false);

  std::stringstream ss;
  ss << "ThermalMasses traces phases in T = [0, Thigh] GeV and prints the "
        "thermal masses.\n"
        "It is called "
        "by\n\n\t./bin/ThermalMasses model input output line"
        "with arguments\n\n\t./bin/ThermalMasses [arguments]\n\nwith the "
        "following arguments, ([*] are required arguments, others "
        "are optional):\n";
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
      arguments.emplace_back("--line=" + std::string(argv[4]));
    }
  }
  return arguments;
}