// Copyright (C) 2024 Lisa Biermann, Margarete Mühlleitner, Rui Santos, João
// Viana SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and
// Jonas Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

/**
 * @file
 * This program calculates Veff on a user-specified field grid
 */

#include <BSMPT/minimizer/Minimizer.h>
#include <BSMPT/minimum_tracer/minimum_tracer.h> // MinimumTracer
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
  int Line{2};
  std::string InputFile, OutputFile;
  double Temperature{-1};
  int npoints{100}, npoints1{0}, npoints2{0}, npoints3{0}, npoints4{0},
      npoints5{0}, npoints6{0};
  double low1{-1}, low2{-1}, low3{-1}, low4{-1}, low5{-1}, low6{-1}, high1{-1},
      high2{-1}, high3{-1}, high4{-1}, high5{-1}, high6{-1};
  bool use_slice_plotter = false;
  std::vector<double> min_start, min_end;
  std::vector<double> point;

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

  std::vector<double> sol, start, solPot;

  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(args.Model, SMConstants);

  std::ifstream infile(args.InputFile);
  if (!infile.good())
  {
    Logger::Write(LoggingLevel::Default, "Input file not found ");
    return EXIT_FAILURE;
  }

  Logger::Write(LoggingLevel::ProgDetailed, "Found file");

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

  double temp = args.Temperature;

  if (args.use_slice_plotter)
  {
    if ((args.min_end.size() == args.min_start.size()) and
        (args.min_end.size() == modelPointer->get_nVEV()))
    {
      for (auto x : modelPointer->addLegendVEV())
        outfile << std::setprecision(16) << x << sep;

      outfile << "Veff(v,T)" << sep << "T" << std::endl;

      Logger::Write(LoggingLevel::ProgDetailed,
                    "Evaluating slice between start minimum at (" +
                        vec_to_string(args.min_start) +
                        ") and end minimum at (" + vec_to_string(args.min_end) +
                        ") with " + std::to_string(args.npoints) +
                        " points at T = " + std::to_string(temp) + " GeV.");

      auto grid_points =
          Create1DimGrid(args.min_start, args.min_end, args.npoints);
      for (auto point : grid_points)
      {
        outfile << point << sep
                << modelPointer->VEff(modelPointer->MinimizeOrderVEV(point),
                                      temp)
                << sep << temp << std::endl;
      }
    }
    else
    {
      Logger::Write(
          LoggingLevel::Default,
          "Given false and/or true minimum has dimensions different from "
          "model VEV dimension.");
      return EXIT_FAILURE;
    }
  }
  else
  {
    for (auto x : modelPointer->addLegendVEV())
      outfile << std::setprecision(16) << x << sep;

    for (auto x : modelPointer->addLegendVEV())
    {
      outfile << x.append("_point") << sep;
    }
    outfile << "Veff(v,T)" << sep << "Veff(point,T)" << sep << "T" << std::endl;

    Logger::Write(LoggingLevel::ProgDetailed,
                  "Temperature is set to T = " + std::to_string(temp));

    std::vector<double> vevStart;
    if (args.point.size() == 0)
    {
      vevStart = std::vector<double>(modelPointer->get_nVEV(), 0);
    }
    else
    {
      if (args.point.size() == modelPointer->get_nVEV())
        vevStart = args.point;
      else
      {
        Logger::Write(LoggingLevel::Default,
                      "Given reference point has dimensions different from "
                      "model VEV dimension.");
        return EXIT_FAILURE;
      }
    }
    Logger::Write(LoggingLevel::ProgDetailed,
                  "Potential grid reference point is: " +
                      vec_to_string(vevStart));

    std::vector<double> vevStartIni = vevStart;

    std::vector<int> npoints = {args.npoints1,
                                args.npoints2,
                                args.npoints3,
                                args.npoints4,
                                args.npoints5,
                                args.npoints6};

    for (std::size_t i = 0; i < vevStartIni.size(); i++)
    {
      if (npoints.at(i) > 0)
      {
        vevStartIni.at(i) =
            0; // reset start vector to zero in directions with npoints > 0
      }
    }

    std::vector<std::vector<double>> res_vec_outer, res_vec_inner_1,
        res_vec_inner_2, res_vec_inner_3, res_vec_inner_4, res_vec_inner_5;

    for (std::size_t i = 0; i < modelPointer->get_nVEV(); i++)
    {
      res_vec_outer =
          Create1DimGrid(vevStartIni, i, args.low1, args.high1, args.npoints1);

      if (modelPointer->get_nVEV() == 1)
      {
        for (std::vector<double> a1 : res_vec_outer)
        {
          outfile << a1 << sep;
          outfile << vevStart << sep;
          outfile << modelPointer->VEff(modelPointer->MinimizeOrderVEV(a1),
                                        temp)
                  << sep;
          outfile << modelPointer->VEff(
                         modelPointer->MinimizeOrderVEV(vevStart), temp)
                  << sep << temp << std::endl;
        }
      }
      else if (i + 1 < modelPointer->get_nVEV())
      {
        for (std::vector<double> a1 : res_vec_outer)
        {
          res_vec_inner_1 =
              Create1DimGrid(a1, i + 1, args.low2, args.high2, args.npoints2);

          if (modelPointer->get_nVEV() == 2)
          {
            for (std::vector<double> a2 : res_vec_inner_1)
            {
              outfile << a2 << sep;
              outfile << vevStart << sep;
              outfile << modelPointer->VEff(modelPointer->MinimizeOrderVEV(a2),
                                            temp)
                      << sep;
              outfile << modelPointer->VEff(
                             modelPointer->MinimizeOrderVEV(vevStart), temp)
                      << sep << temp << std::endl;
            }
          }
          else if (i + 2 < modelPointer->get_nVEV())
          {
            for (std::vector<double> a2 : res_vec_inner_1)
            {
              res_vec_inner_2 = Create1DimGrid(
                  a2, i + 2, args.low3, args.high3, args.npoints3);

              if (modelPointer->get_nVEV() == 3)
              {
                for (std::vector<double> a3 : res_vec_inner_2)
                {
                  outfile << a3 << sep;
                  outfile << vevStart << sep;
                  outfile << modelPointer->VEff(
                                 modelPointer->MinimizeOrderVEV(a3), temp)
                          << sep;
                  outfile << modelPointer->VEff(
                                 modelPointer->MinimizeOrderVEV(vevStart), temp)
                          << sep << temp << std::endl;
                }
              }
              else if (i + 3 < modelPointer->get_nVEV())
              {
                for (std::vector<double> a3 : res_vec_inner_2)
                {
                  res_vec_inner_3 = Create1DimGrid(
                      a3, i + 3, args.low4, args.high4, args.npoints4);

                  if (modelPointer->get_nVEV() == 4)
                  {
                    for (std::vector<double> a4 : res_vec_inner_3)
                    {
                      outfile << a4 << sep;
                      outfile << vevStart << sep;
                      outfile << modelPointer->VEff(
                                     modelPointer->MinimizeOrderVEV(a4), temp)
                              << sep;
                      outfile
                          << modelPointer->VEff(
                                 modelPointer->MinimizeOrderVEV(vevStart), temp)
                          << sep << temp << std::endl;
                    }
                  }
                  else if (i + 4 < modelPointer->get_nVEV())
                  {
                    for (std::vector<double> a4 : res_vec_inner_3)
                    {
                      res_vec_inner_4 = Create1DimGrid(
                          a4, i + 4, args.low5, args.high5, args.npoints5);

                      if (modelPointer->get_nVEV() == 5)
                      {
                        for (std::vector<double> a5 : res_vec_inner_4)
                        {
                          outfile << a5 << sep;
                          outfile << vevStart << sep;
                          outfile
                              << modelPointer->VEff(
                                     modelPointer->MinimizeOrderVEV(a5), temp)
                              << sep;
                          outfile
                              << modelPointer->VEff(
                                     modelPointer->MinimizeOrderVEV(vevStart),
                                     temp)
                              << sep << temp << std::endl;
                        }
                      }
                      else if (i + 5 < modelPointer->get_nVEV())
                      {
                        for (std::vector<double> a5 : res_vec_inner_4)
                        {
                          res_vec_inner_5 = Create1DimGrid(
                              a5, i + 5, args.low6, args.high6, args.npoints6);

                          if (modelPointer->get_nVEV() == 6)
                          {
                            for (std::vector<double> a6 : res_vec_inner_5)
                            {
                              outfile << a6 << sep;
                              outfile << vevStart << sep;
                              outfile << modelPointer->VEff(
                                             modelPointer->MinimizeOrderVEV(a6),
                                             temp)
                                      << sep;
                              outfile << modelPointer->VEff(
                                             modelPointer->MinimizeOrderVEV(
                                                 vevStart),
                                             temp)
                                      << sep << temp << std::endl;
                            }
                          }
                          else if (i + 6 < modelPointer->get_nVEV())
                          {
                            throw std::runtime_error(
                                "Error. More than 6-dim grids "
                                "are not implemented.");
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
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
  if (Model == ModelID::ModelIDs::NotSet)
  {
    Logger::Write(LoggingLevel::Default,
                  "Your Model parameter does not match with the "
                  "implemented Models.");
    ShowInputError();
    return false;
  }
  if (Temperature < 0)
  {
    if (Temperature == -1)
      Logger::Write(LoggingLevel::Default, "No input for temperature.");
    else
      Logger::Write(LoggingLevel::Default, "Invalid input for temperature.");
    return false;
  }
  if (use_slice_plotter and ((min_start.size() == 0 or min_end.size() == 0)))
  {
    Logger::Write(
        LoggingLevel::Default,
        "Slice plotter mode chosen but no input for false and/or true vacuum.");
    return false;
  }
  if ((npoints1 < 0) or (npoints2 < 0) or (npoints3 < 0) or (npoints4 < 0) or
      (npoints5 < 0) or (npoints6 < 0) or (npoints < 0))
  {
    Logger::Write(LoggingLevel::Default, "Invalid grid size requested.");
    return false;
  }
  if ((low1 > high1) or (low2 > high2) or (low3 > high3) or (low4 > high4) or
      (low5 > high5) or (low6 > high6))
  {
    Logger::Write(LoggingLevel::Default,
                  "Invalid field value boundaries requested.");
    return false;
  }
  if ((npoints1 > 0 and (low1 == -1 or high1 == -1)) or
      (npoints2 > 0 and (low2 == -1 or high2 == -1)) or
      (npoints3 > 0 and (low3 == -1 or high3 == -1)) or
      (npoints4 > 0 and (low4 == -1 or high4 == -1)) or
      (npoints5 > 0 and (low5 == -1 or high5 == -1)) or
      (npoints6 > 0 and (low6 == -1 or high6 == -1)))
  {
    Logger::Write(LoggingLevel::Default,
                  "Grid size in VEV direction non-zero, but no low/high field "
                  "values defined.");
    return false;
  }
  return true;
}

CLIOptions::CLIOptions(const BSMPT::parser &argparser)
{
  std::stringstream ss;
  argparser.check_required_parameters();

  // required arguments
  Model       = BSMPT::ModelID::getModel(argparser.get_value("model"));
  InputFile   = argparser.get_value("input");
  OutputFile  = argparser.get_value("output");
  Line        = argparser.get_value<int>("line");
  Temperature = argparser.get_value<double>("temperature");

  // optional arguments

  try
  {
    auto vec_str = split(argparser.get_value("point"), ',');
    for (std::size_t i = 0; i < vec_str.size(); i++)
    {
      point.push_back(std::stod(vec_str.at(i)));
    }
  }
  catch (BSMPT::parserException &)
  {
    ss << "--point not set\n";
  }

  try
  {
    npoints1 = argparser.get_value<int>("npoints1");
  }
  catch (BSMPT::parserException &)
  {
    ss << "--npoints1 not set, using default value: " << npoints1 << "\n";
  }
  try
  {
    npoints2 = argparser.get_value<int>("npoints2");
  }
  catch (BSMPT::parserException &)
  {
    ss << "--npoints2 not set, using default value: " << npoints2 << "\n";
  }
  try
  {
    npoints3 = argparser.get_value<int>("npoints3");
  }
  catch (BSMPT::parserException &)
  {
    ss << "--npoints3 not set, using default value: " << npoints3 << "\n";
  }
  try
  {
    npoints4 = argparser.get_value<int>("npoints4");
  }
  catch (BSMPT::parserException &)
  {
    ss << "--npoints4 not set, using default value: " << npoints4 << "\n";
  }
  try
  {
    npoints5 = argparser.get_value<int>("npoints5");
  }
  catch (BSMPT::parserException &)
  {
    ss << "--npoints5 not set, using default value: " << npoints5 << "\n";
  }
  try
  {
    npoints6 = argparser.get_value<int>("npoints6");
  }
  catch (BSMPT::parserException &)
  {
    ss << "--npoints6 not set, using default value: " << npoints6 << "\n";
  }

  try
  {
    low1 = argparser.get_value<double>("low1");
  }
  catch (BSMPT::parserException &)
  {
    ss << "--low1 not set, using default value: " << low1 << "\n";
  }
  try
  {
    low2 = argparser.get_value<double>("low2");
  }
  catch (BSMPT::parserException &)
  {
    ss << "--low2 not set, using default value: " << low2 << "\n";
  }
  try
  {
    low3 = argparser.get_value<double>("low3");
  }
  catch (BSMPT::parserException &)
  {
    ss << "--low3 not set, using default value: " << low3 << "\n";
  }
  try
  {
    low4 = argparser.get_value<double>("low4");
  }
  catch (BSMPT::parserException &)
  {
    ss << "--low4 not set, using default value: " << low4 << "\n";
  }
  try
  {
    low5 = argparser.get_value<double>("low5");
  }
  catch (BSMPT::parserException &)
  {
    ss << "--low5 not set, using default value: " << low5 << "\n";
  }
  try
  {
    low6 = argparser.get_value<double>("low6");
  }
  catch (BSMPT::parserException &)
  {
    ss << "--low6 not set, using default value: " << low6 << "\n";
  }

  try
  {
    high1 = argparser.get_value<double>("high1");
  }
  catch (BSMPT::parserException &)
  {
    ss << "--high1 not set, using default value: " << high1 << "\n";
  }
  try
  {
    high2 = argparser.get_value<double>("high2");
  }
  catch (BSMPT::parserException &)
  {
    ss << "--high2 not set, using default value: " << high2 << "\n";
  }
  try
  {
    high3 = argparser.get_value<double>("high3");
  }
  catch (BSMPT::parserException &)
  {
    ss << "--high3 not set, using default value: " << high3 << "\n";
  }
  try
  {
    high4 = argparser.get_value<double>("high4");
  }
  catch (BSMPT::parserException &)
  {
    ss << "--high4 not set, using default value: " << high4 << "\n";
  }
  try
  {
    high5 = argparser.get_value<double>("high5");
  }
  catch (BSMPT::parserException &)
  {
    ss << "--high5 not set, using default value: " << high5 << "\n";
  }
  try
  {
    high6 = argparser.get_value<double>("high6");
  }
  catch (BSMPT::parserException &)
  {
    ss << "--high6 not set, using default value: " << high6 << "\n";
  }

  try
  {
    use_slice_plotter = (argparser.get_value("slice") == "true");
  }
  catch (BSMPT::parserException &)
  {
    ss << "--slice not set, using default value: false \n";
  }

  try
  {
    auto vec_str = split(argparser.get_value("min_start"), ',');
    for (std::size_t i = 0; i < vec_str.size(); i++)
    {
      min_start.push_back(std::stod(vec_str.at(i)));
    }
  }
  catch (BSMPT::parserException &)
  {
    ss << "--min_start not set\n";
  }

  try
  {
    auto vec_str = split(argparser.get_value("min_end"), ',');
    for (std::size_t i = 0; i < vec_str.size(); i++)
    {
      min_end.push_back(std::stod(vec_str.at(i)));
    }
  }
  catch (BSMPT::parserException &)
  {
    ss << "--min_end not set\n";
  }

  try
  {
    npoints = argparser.get_value<int>("npoints");
  }
  catch (BSMPT::parserException &)
  {
    ss << "--npoints not set, using default value: " << npoints << "\n";
  }
}

BSMPT::parser prepare_parser()
{
  BSMPT::parser argparser(true);
  argparser.add_argument("help", "shows this menu", false);
  argparser.add_argument("model", "[*] model name", true);
  argparser.add_argument("input", "[*] input file (in tsv format)", true);
  argparser.add_argument("output", "[*] output file (in tsv format)", true);
  argparser.add_argument("line", "[*] line number of line in input file", true);
  argparser.add_subtext("    (expects line 1 to be a legend)");
  argparser.add_argument("temperature", "[*] temperature [GeV]", true);
  argparser.add_argument("point", "grid reference point", "0,..,0", false);
  argparser.add_argument_only_display(
      "npointsi", "number of points in direction i", "0");
  argparser.add_subtext("(with i = [1,..,6])");
  argparser.add_argument("npoints1", false);
  argparser.add_argument("npoints2", false);
  argparser.add_argument("npoints3", false);
  argparser.add_argument("npoints4", false);
  argparser.add_argument("npoints5", false);
  argparser.add_argument("npoints6", false);
  argparser.add_argument_only_display(
      "lowi", "lowest field value in direction i", "0");
  argparser.add_subtext("[* if npointsi > 0] (with i = [1,..,6])");
  argparser.add_argument("low1", false);
  argparser.add_argument("low2", false);
  argparser.add_argument("low3", false);
  argparser.add_argument("low4", false);
  argparser.add_argument("low5", false);
  argparser.add_argument("low6", false);
  argparser.add_argument_only_display(
      "highi", "highest field value in direction i", "0");
  argparser.add_subtext("[* if npointsi > 0] (with i = [1,..,6])");
  argparser.add_argument("high1", false);
  argparser.add_argument("high2", false);
  argparser.add_argument("high3", false);
  argparser.add_argument("high4", false);
  argparser.add_argument("high5", false);
  argparser.add_argument("high6", false);
  argparser.add_argument("slice", "enable slice mode", "false", false);
  argparser.add_argument("min_start", "[* in slice mode] start minimum", false);
  argparser.add_argument("min_end", "[* in slice mode] end minimum", false);
  argparser.add_argument("npoints", "grid size in slice mode", "100", false);
  argparser.add_argument(
      "json", "use a json file instead of cli parameters", "", false);

  std::stringstream ss;
  ss << "PotPlotter calculates the effective potential on a user-specified "
        "field grid\nit is "
        "called "
        "by\n\n\t./bin/PotPlotter [arguments]\n\nwith the "
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
  return arguments;
}