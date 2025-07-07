// Copyright (C) 2024 Lisa Biermann, Margarete Mühlleitner, Rui Santos, João
// Viana SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and
// Jonas Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

/**
 * @file
 * This program calculates characteristic temperatures for
 * phase transitions
 *
 */

#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/transition_tracer/transition_tracer.h>
#include <BSMPT/utility/Logger.h>
#include <BSMPT/utility/parser.h>
#include <BSMPT/utility/utility.h>
#include <Eigen/Dense>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdlib.h> // for atoi, EXIT_FAILURE
#include <string>   // for string, operator<<
#include <utility>  // for pair
#include <vector>   // for vector

using namespace std;
using namespace BSMPT;

struct CLIOptions
{
  BSMPT::ModelID::ModelIDs Model{ModelID::ModelIDs::NotSet};
  int firstline{0}, lastline{0};
  double templow{0}, temphigh{300};
  double UserDefined_vwall = 0.95;
  int MaxPathIntegrations  = 7;
  std::string inputfile, outputfile;
  bool UseGSL{Minimizer::UseGSLDefault};
  bool UseCMAES{Minimizer::UseLibCMAESDefault};
  bool UseNLopt{Minimizer::UseNLoptDefault};
  int WhichMinimizer{Minimizer::WhichMinimizerDefault};
  bool UseMultithreading{false};
  int UseMultiStepPTMode{-1};
  int CheckEWSymmetryRestoration{1};
  double perc_prbl{.71};
  double compl_prbl{.01};
  int num_check_pts{10};
  int CheckNLOStability{1};

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
  int linecounter   = 1;
  std::size_t count = 0;
  int num_points    = args.lastline - args.firstline + 1;

  // output contents storage
  std::vector<std::stringstream> output_contents;
  output_contents.resize(num_points); // reserve one row per point
  std::vector<std::string> transition_history;
  std::vector<std::string> legend;

  while (getline(infile, linestr))
  {
    if (linecounter == 1) linestr_store = linestr;

    if (linecounter > args.lastline) break;

    if (linecounter >= args.firstline and linecounter <= args.lastline)
    {
      count += 1; // keep track at which point we are
      output_contents.at(count - 1).precision(
          std::numeric_limits<double>::max_digits10);

      Logger::Write(LoggingLevel::ProgDetailed,
                    "Currently at line " + std::to_string(linecounter));

      modelPointer->setUseIndexCol(linestr_store);

      std::pair<std::vector<double>, std::vector<double>> parameters =
          modelPointer->initModel(linestr);

      if (args.firstline == args.lastline)
      {
        modelPointer->write();
      }

      auto start = std::chrono::high_resolution_clock::now();

      user_input input{modelPointer,
                       args.templow,
                       args.temphigh,
                       args.UserDefined_vwall,
                       args.perc_prbl,
                       args.compl_prbl,
                       0.1,
                       args.MaxPathIntegrations,
                       args.UseMultiStepPTMode,
                       args.num_check_pts,
                       args.CheckEWSymmetryRestoration,
                       args.CheckNLOStability,
                       args.WhichMinimizer,
                       args.UseMultithreading,
                       false,
                       TransitionTemperature::Percolation,
                       1};

      TransitionTracer trans(input);

      auto time = std::chrono::duration_cast<std::chrono::milliseconds>(
                      std::chrono::high_resolution_clock::now() - start)
                      .count() /
                  1000.;

      BSMPT::Logger::Write(BSMPT::LoggingLevel::ProgDetailed,
                           "\nTook\t" + std::to_string(time) + " seconds.\n");

      auto output = trans.output_store;

      output_contents.at(count - 1)
          << linestr << sep << parameters.second << sep
          << output.status.status_nlo_stability << sep
          << output.status.status_ewsr << sep << output.status.status_tracing
          << sep << output.status.status_coex_pairs << sep << time << sep;

      if ((output.status.status_tracing == StatusTracing::Success) &&
          (output.status.status_coex_pairs == StatusCoexPair::Success))
      {
        for (std::size_t i = 0; i < trans.output_store.num_coex_phase_pairs;
             i++)
        {
          output_contents.at(count - 1)
              << output.status.status_crit.at(i) << sep
              << output.vec_trans_data.at(i).crit_temp.value_or(EmptyValue)
              << sep << output.vec_trans_data.at(i).crit_false_vev << sep
              << output.vec_trans_data.at(i).crit_true_vev << sep
              << output.status.status_bounce_sol.at(i) << sep
              << output.status.status_nucl_approx.at(i) << sep
              << output.vec_trans_data.at(i).nucl_approx_temp.value_or(
                     EmptyValue)
              << sep << output.vec_trans_data.at(i).nucl_approx_false_vev << sep
              << output.vec_trans_data.at(i).nucl_approx_true_vev << sep
              << output.status.status_nucl.at(i) << sep
              << output.vec_trans_data.at(i).nucl_temp.value_or(EmptyValue)
              << sep << output.vec_trans_data.at(i).nucl_false_vev << sep
              << output.vec_trans_data.at(i).nucl_true_vev << sep
              << output.status.status_perc.at(i) << sep
              << output.vec_trans_data.at(i).perc_temp.value_or(EmptyValue)
              << sep << output.vec_trans_data.at(i).perc_false_vev << sep
              << output.vec_trans_data.at(i).perc_true_vev << sep
              << output.status.status_compl.at(i) << sep
              << output.vec_trans_data.at(i).compl_temp.value_or(EmptyValue)
              << sep << output.vec_trans_data.at(i).compl_false_vev << sep
              << output.vec_trans_data.at(i).compl_true_vev << sep;
        }
      }

      transition_history.push_back(output.transition_history);

      if (legend.size() < output.legend.size())
      {
        legend = output.legend; // update legend
      }

      // write to output file
      std::ofstream outfile(args.outputfile);
      if (!outfile.good())
      {
        Logger::Write(LoggingLevel::Default,
                      "Can not create file " + args.outputfile);
        return EXIT_FAILURE;
      }

      int tab_count_legend = 1;
      std::stringstream full_legend;
      full_legend << linestr_store << sep << modelPointer->addLegendCT() << sep
                  << legend;
      outfile << full_legend.str() << std::endl;

      // get length of legend and of contents
      std::string str_legend = full_legend.str();
      for (auto &el : str_legend)
      {
        if (el == '\t')
        {
          tab_count_legend += 1;
        }
      }

      std::vector<int> tab_count(count, 1);
      for (std::size_t i = 0; i < count; i++)
      {
        std::string line = output_contents.at(i).str();
        for (auto &el : line)
        {
          if (el == '\t')
          {
            tab_count.at(i) += 1;
          }
        }
      }

      // fill up previous rows to match multi-line
      for (std::size_t i = 0; i < count; i++)
      {
        outfile << output_contents.at(i).str();

        int diff = int(tab_count_legend - tab_count.at(i));

        while (diff > 0)
        {
          outfile << "nan" << sep;
          diff--;
        }
        outfile << transition_history.at(i) << std::endl;
      }

      outfile.close();
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
    throw std::runtime_error("You set --UseGSL=true but GSL was not "
                             "found during compilation.");
  }
  if (UseCMAES and not Minimizer::UseLibCMAESDefault)
  {
    throw std::runtime_error("You set --UseCMAES=true but CMAES was not "
                             "found during compilation.");
  }
  if (UseNLopt and not Minimizer::UseNLoptDefault)
  {
    throw std::runtime_error("You set --UseNLopt=true but NLopt was not "
                             "found during compilation.");
  }
  if (WhichMinimizer == 0)
  {
    throw std::runtime_error(
        "You disabled all minimizers. You need at least one.");
  }

  if (Model == ModelID::ModelIDs::NotSet)
  {
    Logger::Write(LoggingLevel::Default,
                  "Your Model parameter does not match with the "
                  "implemented Models.");
    ShowInputError();
    return false;
  }
  if (firstline == 0 or lastline == 0)
  {
    Logger::Write(LoggingLevel::Default, "firstline or lastline not set.");
    return false;
  }
  if (firstline < 0 or lastline < 0)
  {
    Logger::Write(LoggingLevel::Default,
                  "Invalid input for first- or lastline.");
    return false;
  }
  if (firstline > lastline)
  {
    Logger::Write(LoggingLevel::Default, "lastline is smaller then firstline.");
    return false;
  }
  if (templow >= temphigh)
  {
    Logger::Write(LoggingLevel::Default,
                  "Invalid temperature choice. Thigh has to be > 0 GeV.");
    return false;
  }
  if ((UserDefined_vwall != -1 and UserDefined_vwall != -2) and
      (UserDefined_vwall <= 0 or UserDefined_vwall > 1))
  {
    Logger::Write(LoggingLevel::Default, "Invalid choice for vwall.");
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
  if (CheckNLOStability < 0 or CheckNLOStability > 1)
  {
    Logger::Write(LoggingLevel::Default,
                  "Invalid choice for CheckNLOStability.");
    return false;
  }
  if (perc_prbl < 0 or perc_prbl > 1)
  {
    Logger::Write(LoggingLevel::Default,
                  "Invalid false vacuum fraction for determination of "
                  "percolation temperature given.");
    return false;
  }
  if (compl_prbl < 0 or compl_prbl > 1)
  {
    Logger::Write(LoggingLevel::Default,
                  "Invalid false vacuum fraction for determination of "
                  "completion temperature given.");
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
  firstline  = argparser.get_value<int>("firstline");
  lastline   = argparser.get_value<int>("lastline");

  // optional arguments
  try
  {
    temphigh = argparser.get_value<double>("thigh");
  }
  catch (BSMPT::parserException &)
  {
    ss << "--thigh not set, using default value: " << temphigh << "\n";
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

  try
  {
    UserDefined_vwall = argparser.get_value<double>("vwall");
  }
  catch (BSMPT::parserException &)
  {
    ss << "--vwall not set, using default value: " << UserDefined_vwall << "\n";
  }

  try
  {
    perc_prbl = argparser.get_value<double>("perc_prbl");
  }
  catch (BSMPT::parserException &)
  {
    ss << "--perc_prbl not set, using default value: " << perc_prbl << "\n";
  }

  try
  {
    compl_prbl = argparser.get_value<double>("compl_prbl");
  }
  catch (BSMPT::parserException &)
  {
    ss << "--compl_prbl not set, using default value: " << compl_prbl << "\n";
  }

  // CheckNLOStability
  try
  {
    auto nlo_string = argparser.get_value("checknlo");
    if (nlo_string == "on")
    {
      CheckNLOStability = 1;
    }
    else if (nlo_string == "off")
    {
      CheckNLOStability = 0;
    }
  }
  catch (BSMPT::parserException &)
  {
    ss << "--checknlo not set, using default value: on\n";
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
    else if (ewsr_string == "keep_bfb")
    {
      CheckEWSymmetryRestoration = 2;
    }
    else if (ewsr_string == "keep_ewsr")
    {
      CheckEWSymmetryRestoration = 3;
    }
  }
  catch (BSMPT::parserException &)
  {
    ss << "--checkewsr not set, using default value: on\n";
  }

  try
  {
    MaxPathIntegrations = argparser.get_value<int>("maxpathintegrations");
  }
  catch (BSMPT::parserException &)
  {
    ss << "--maxpathintegrations not set, using default value: "
       << MaxPathIntegrations << "\n";
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
      "firstline", "[*] line number of first line in input file", true);
  argparser.add_subtext("    (expects line 1 to be a legend)");
  argparser.add_argument(
      "lastline", "[*] line number of last line in input file", true);
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
      "vwall", "wall velocity: >0 user defined", "0.95", false);
  argparser.add_subtext("-1: approximation");
  argparser.add_subtext("-2: upper bound");
  argparser.add_argument(
      "perc_prbl", "false vacuum fraction for percolation", "0.71", false);
  argparser.add_argument(
      "compl_prbl", "false vacuum fraction for completion", "0.01", false);
  argparser.add_argument("checknlo", "check for NLO stability", "on", false);
  argparser.add_subtext("on: only keep NLO stable points");
  argparser.add_subtext("off: check disabled");
  argparser.add_argument(
      "checkewsr", "check for EWSR at high temperature", "on", false);
  argparser.add_subtext("on: perform check and add info");
  argparser.add_subtext("keep_bfb: only keep BFB points");
  argparser.add_subtext("keep_ewsr: only keep EWSR points");
  argparser.add_subtext("off: check disabled");

  argparser.add_argument("maxpathintegrations",
                         "number of solutions of 1D equation =",
                         "7",
                         false);
  argparser.add_subtext("number of path deformations + 1");

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
  ss << "CalcTemps calculates characteristic temperatures for phase "
        "transitions\nit is called "
        "by\n\n\t./bin/CalcTemps model input output firstline "
        "lastline\n\nor "
        "with arguments\n\n\t./bin/CalcTemps [arguments]\n\nwith the "
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
      arguments.emplace_back("--firstline=" + std::string(argv[4]));
    }
    if (argc >= 6)
    {
      arguments.emplace_back("--lastline=" + std::string(argv[5]));
    }
  }
  return arguments;
}