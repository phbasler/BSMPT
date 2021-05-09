// Copyright (C) 2020  Philipp Basler, Margarete Mühlleitner and Jonas
// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

/**
 * @file
 * Calculates the Counterterms for mulitple points and adds them at the end of
 * the line. The sequence is the same as used in Potential::set_CT_Pot_Par
 */

#include <BSMPT/models/ClassPotentialOrigin.h> // for Class_Potential_Origin
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/utility/Logger.h>
#include <BSMPT/utility/utility.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>   // for unique_ptr
#include <stdlib.h> // for atoi, EXIT_FAILURE
#include <string>   // for operator<<, getline
#include <utility>  // for pair
#include <vector>   // for vector
using namespace std;
using namespace BSMPT;

struct CLIOptions
{
  BSMPT::ModelID::ModelIDs Model{ModelID::ModelIDs::NotSet};
  int FirstLine{}, LastLine{};
  std::string InputFile, OutputFile;
  bool TerminalOutput{false};

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

  std::unique_ptr<Class_Potential_Origin> modelPointer =
      ModelID::FChoose(args.Model);

  while (getline(infile, linestr))
  {
    if (linecounter > args.LastLine) break;
    if (linecounter == 1)
    {
      modelPointer->setUseIndexCol(linestr);
      outfile << linestr;
      outfile << sep << modelPointer->addLegendCT();
      outfile << std::endl;
    }

    if (linecounter >= args.FirstLine and linecounter <= args.LastLine and
        linecounter != 1)
    {
      std::pair<std::vector<double>, std::vector<double>> parameters =
          modelPointer->initModel(linestr);
      outfile << linestr << sep << parameters.second << std::endl;
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
    ss << "CalcCT calculates the counterterm parameters for the given "
          "paramter points"
       << std::endl
       << "It is called either by " << std::endl
       << "./CalcCT model input output FirstLine LastLine" << std::endl
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
       << "The first line in the input file to calculate the EWPT. "
          "Expects line 1 to be a legend."
       << std::endl
       << std::setw(SizeOfFirstColumn) << std::left << "--LastLine="
       << "The last line in the input file to calculate the EWPT." << std::endl
       << std::setw(SizeOfFirstColumn) << std::left << "--TerminalOutput="
       << "y/n Turns on additional information in the terminal during "
          "the calculation."
       << std::endl;
    Logger::Write(LoggingLevel::Default, ss.str());
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

  std::string prefix{"--"};
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
      else if (StringStartsWith(el, "--firstline="))
      {
        FirstLine = std::stoi(el.substr(std::string("--firstline=").size()));
      }
      else if (StringStartsWith(el, "--lastline="))
      {
        LastLine = std::stoi(el.substr(std::string("--lastline=").size()));
      }
      else if (StringStartsWith(el, "--terminaloutput="))
      {
        TerminalOutput =
            el.substr(std::string("--terminaloutput=").size()) == "y";
      }
    }
  }
  else
  {
    Model      = ModelID::getModel(args.at(0));
    InputFile  = args.at(1);
    OutputFile = args.at(2);
    FirstLine  = std::stoi(args.at(3));
    LastLine   = std::stoi(args.at(4));
    if (argc == 7)
    {
      TerminalOutput = ("y" == std::string(argv[6]));
    }
  }
}

bool CLIOptions::good() const
{
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
