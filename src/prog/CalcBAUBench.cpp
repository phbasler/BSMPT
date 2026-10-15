// Copyright (C) 2024 Lisa Biermann, Margarete Mühlleitner, Rui Santos, João
// Viana SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and
// Jonas Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

/**
 * @file
 * This program calculates properties of graviational waves sourced by phase
 * transitions
 *
 */

#include <BSMPT/baryo_fhck/BenchmarkModel.h>
#include <BSMPT/baryo_fhck/TransportEquations.h>
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/transition_tracer/transition_tracer.h>
#include <BSMPT/utility/Logger.h>
#include <BSMPT/utility/parser.h>
#include <BSMPT/utility/utility.h>
#include <Eigen/Dense>
#include <chrono>
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
  std::string inputfile, outputfile;
  std::vector<size_t> FHCKMoments = {2};
  BSMPT::Baryo::FHCK::TruncationScheme truncationscheme =
      BSMPT::Baryo::FHCK::TruncationScheme::MinusVw;
  double truncationR = 0;
  CLIOptions(const BSMPT::parser &argparser);
  bool good() const;
};

BSMPT::parser prepare_parser();

std::vector<std::string> convert_input(int argc, char *argv[]);

int main(int argc, char *argv[])
try
{
  auto argparser = prepare_parser();
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

  std::string linestr, linestr_store;
  int linecounter   = 1;
  std::size_t count = 0;
  int num_points    = args.lastline - args.firstline + 1;

  // output contents storage
  std::vector<std::stringstream> output_contents;
  output_contents.resize(num_points); // reserve one row per point

  std::vector<std::string> legend;

  while (getline(infile, linestr))
  {
    if (linecounter == 1)
    {
      linestr_store = linestr;
    }

    if (linecounter > args.lastline) break;

    if (linecounter >= args.firstline and linecounter <= args.lastline)
    {
      count += 1; // keep track at which point we are
      output_contents.at(count - 1).precision(
          std::numeric_limits<double>::max_digits10);

      Logger::Write(LoggingLevel::ProgDetailed,
                    "Currently at line " + std::to_string(linecounter));

      user_input input{};

      // fallback -> vw = 0.95
      /* const double vwall_with_fallback =
          (args.UserDefined_vwall < 0 and not input.only_crit)
              ? output.vec_gw_data.at(i).vwall.value()
              : args.UserDefined_vwall; */
      output_contents.at(count - 1) << linestr << sep;
      std::stringstream ss(linestr);
      double vn, wn, Tn, LAM, Lw, Ls, vw;
      ss >> vn >> wn >> Tn >> LAM >> Lw >> Ls >> vw;

      std::shared_ptr<BSMPT::Baryo::FHCK::BenchmarkModel> transportmodel =
          std::make_shared<BSMPT::Baryo::FHCK::BenchmarkModel>(
              vn,
              wn,
              Tn,
              LAM,
              Lw,
              Ls,
              vw,
              args.truncationscheme,
              args.truncationR);

      BSMPT::Baryo::FHCK::TransportEquations transportequation(
          transportmodel, Tn, args.FHCKMoments);
      transportequation.SolveTransportEquation();

      std::stringstream full_legend;
      full_legend << linestr_store;

      for (const auto &eta : transportequation.BAUeta)
      {
        output_contents.at(count - 1) << eta.value_or(EmptyValue) << sep;
        full_legend << sep << "eta("
                    << args.FHCKMoments.at(&eta -
                                           &transportequation.BAUeta.at(0))
                    << ")";
      }

      // write to output file
      std::ofstream outfile(args.outputfile);
      if (!outfile.good())
      {
        Logger::Write(LoggingLevel::Default,
                      "Can not create file " + args.outputfile);
        return EXIT_FAILURE;
      }
      outfile << full_legend.str() << std::endl;

      // fill up previous rows to match multi-line
      for (std::size_t i = 0; i < count; i++)
        outfile << output_contents.at(i).str() << std::endl;

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

CLIOptions::CLIOptions(const BSMPT::parser &argparser)
{
  std::stringstream ss;
  argparser.check_required_parameters();
  inputfile  = argparser.get_value("input");
  outputfile = argparser.get_value("output");
  firstline  = argparser.get_value<int>("firstline");
  lastline   = argparser.get_value<int>("lastline");

  try
  {
    auto truncationschemestring = argparser.get_value("truncationscheme");

    if (truncationschemestring == "minusvw")
    {
      truncationscheme = BSMPT::Baryo::FHCK::TruncationScheme::MinusVw;
    }
    else if (truncationschemestring == "variance")
    {
      truncationscheme = BSMPT::Baryo::FHCK::TruncationScheme::Variance;
    }
    else
    {
      try
      {
        truncationR      = std::stod(truncationschemestring);
        truncationscheme = BSMPT::Baryo::FHCK::TruncationScheme::Const;
      }
      catch (BSMPT::parserException &)
      {
        ss << "--truncationscheme set with invalid option: '"
           << truncationschemestring << "'. Defaulting back to 'minusvw'.\n";
        truncationscheme = BSMPT::Baryo::FHCK::TruncationScheme::MinusVw;
      }
    }
  }
  catch (BSMPT::parserException &)
  {
    truncationscheme = BSMPT::Baryo::FHCK::TruncationScheme::MinusVw;
    ss << "--truncationscheme not set, using default value: 'minusvw' \n";
  }

  try
  {
    auto vec_str = split(argparser.get_value("moments"), ',');
    FHCKMoments.clear();
    for (std::size_t i = 0; i < vec_str.size(); i++)
    {
      size_t moment = std::stoi(vec_str.at(i));
      if ((moment - 2) % 4 == 0)
        FHCKMoments.push_back(moment); // check if 2 + 4n
    }
    // sort and remove duplicates
    std::sort(FHCKMoments.begin(), FHCKMoments.end());
    FHCKMoments.erase(unique(FHCKMoments.begin(), FHCKMoments.end()),
                      FHCKMoments.end());
  }
  catch (BSMPT::parserException &)
  {
    ss << "--moments not set, using default value: " << FHCKMoments << "\n";
  }
}

bool CLIOptions::good() const
{
  if (firstline == 0 or lastline == 0)
  {
    Logger::Write(LoggingLevel::Default, "firstline or lastline not set.");
    return false;
  }
  if (firstline > lastline)
  {
    Logger::Write(LoggingLevel::Default, "firstline is smaller then lastline.");
    return false;
  }
  return true;
}

BSMPT::parser prepare_parser()
{
  BSMPT::parser argparser(true);
  argparser.add_argument("help", "shows this menu", false);
  argparser.add_argument("input", "[*] input file (in tsv format)", true);
  argparser.add_argument("output", "[*] output file (in tsv format)", true);
  argparser.add_argument(
      "firstline", "[*] line number of first line in input file", true);
  argparser.add_subtext("    (expects line 1 to be a legend)");
  argparser.add_argument(
      "lastline", "[*] line number of last line in input file", true);
  argparser.add_argument(
      "truncationscheme", "truncation scheme to be used", "minusvw", false);
  argparser.add_subtext("minusvw: R = -vw");
  argparser.add_subtext("<number>: R = number");
  argparser.add_subtext("variance: variance truncation");
  argparser.add_argument(
      "moments", "moments to solve the transport equations", "2", false);

  std::stringstream ss;
  ss << "CalcBAUBench calculates the baryon assymetry of the Universe for the "
        "benchmark model in [2407.13639] \nit is "
        "called "
        "by\n\n\t./bin/CalcBAUBench input output firstline "
        "lastline\n\nor "
        "with arguments\n\n\t./bin/CalcBAUBench[arguments]\n\nwith the "
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
      arguments.emplace_back("--input=" + std::string(argv[1]));
    }
    if (argc >= 3)
    {
      arguments.emplace_back("--output=" + std::string(argv[2]));
    }
    if (argc >= 4)
    {
      arguments.emplace_back("--firstline=" + std::string(argv[3]));
    }
    if (argc >= 5)
    {
      arguments.emplace_back("--lastline=" + std::string(argv[4]));
    }
  }
  return arguments;
}
