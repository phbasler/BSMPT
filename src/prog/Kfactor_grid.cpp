// Copyright (C) 2020  Philipp Basler, Margarete Mühlleitner and Jonas
// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

/**
 * @file
 * Calculates the electroweak phase transition for a given Inputfile for a given
 * subset of lines in the file and adds it at the end of the line in the format
 * T_c v_c all single vevs. One parameter point per line.
 *
 */

#include <BSMPT/Kfactors/Kfactors.h>
#include <BSMPT/models/IncludeAllModels.h>
#include <algorithm> // for max, copy
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory> // for allocator_traits<>::value...
#include <sstream>
#include <stdlib.h> // for std::size_t, atoi, EXIT_FAILURE
#include <string>   // for string, operator==, opera...
#include <vector>   // for vector
using namespace std;
using namespace BSMPT;

void CreateKtildeInterpolationData();
void CreateKfunctionsgrid();

void CreateKtildeInterpolationData()
{
  std::vector<double> msquaredTsquared, KtildeNormBoson, KtildeNormFermion;
  for (int i = 0; i <= 400; i++)
  {
    msquaredTsquared.push_back(i);
    KtildeNormBoson.push_back(Kfactors::Ktilde_normalization(1, -1, i));
    KtildeNormFermion.push_back(Kfactors::Ktilde_normalization(1, 1, i));
  }

  std::ofstream outputNormalisation("KtildeInterpolation.h");

  outputNormalisation
      << "/*"
      << "\n"
      << " * KtildeInterpolation.h"
      << "\n"
      << " *"
      << "\n"
      << " *  Copyright (C) 2020  Philipp Basler, Margarete Mühlleitner and "
         "Jonas Müller"
      << "\n"
      << ""
      << "\n"
      << "      This program is free software: you can redistribute it and/or "
         "modify"
      << "\n"
      << "      it under the terms of the GNU General Public License as "
         "published by"
      << "\n"
      << "      the Free Software Foundation, either version 3 of the License, "
         "or"
      << "\n"
      << "      (at your option) any later version."
      << "\n"
      << ""
      << "\n"
      << "      This program is distributed in the hope that it will be useful,"
      << "\n"
      << "      but WITHOUT ANY WARRANTY; without even the implied warranty of"
      << "\n"
      << "      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the"
      << "\n"
      << "      GNU General Public License for more details."
      << "\n"
      << ""
      << "\n"
      << "      You should have received a copy of the GNU General Public "
         "License"
      << "\n"
      << "      along with this program.  If not, see "
         "<http://www.gnu.org/licenses/>."
      << "\n"
      << " */"
      << "\n";

  outputNormalisation
      << "#ifndef SRC_BARYO_CALCULATION_KFACTORS_GRID_KtildeInterpolation_H_ "
      << "\n"
      << "#define SRC_BARYO_CALCULATION_KFACTORS_GRID_KtildeInterpolation_H_"
      << "\n";

  outputNormalisation << std::endl;

  outputNormalisation
      << "/**"
      << "\n"
      << "* @file"
      << "\n"
      << "* Data points for the interpolation of the normalisation of Ktile"
      << "\n"
      << "*/" << std::endl;

  outputNormalisation << "namespace BSMPT { \n"
                      << "namespace Kfactors { \n"
                      << "namespace Data{ \n"
                      << std::endl;

  outputNormalisation << "const int KtildeInterpolationSize = "
                      << msquaredTsquared.size() << ";" << std::endl;

  outputNormalisation
      << "const std::vector<double> msquaredTsquared{"; //<<
                                                        //msquaredTsquared.size()
                                                        //<< "] = {";
  for (std::size_t i = 0; i < msquaredTsquared.size(); i++)
  {
    outputNormalisation << msquaredTsquared.at(i);
    if (i < msquaredTsquared.size() - 1) outputNormalisation << ",";
  }
  outputNormalisation << "};" << std::endl;

  outputNormalisation
      << "const std::vector<double> KtildeNormBoson_grid{"; // <<
                                                            // KtildeNormBoson.size()
                                                            // << "] = {";
  for (std::size_t i = 0; i < KtildeNormBoson.size(); i++)
  {
    outputNormalisation << KtildeNormBoson.at(i);
    if (i < KtildeNormBoson.size() - 1) outputNormalisation << ",";
  }
  outputNormalisation << "};" << std::endl;

  outputNormalisation
      << "const std::vector<double> KtildeNormFermion_grid{"; //<<
                                                              //KtildeNormFermion.size()
                                                              //<< "] = {";
  for (std::size_t i = 0; i < KtildeNormFermion.size(); i++)
  {
    outputNormalisation << KtildeNormFermion.at(i);
    if (i < KtildeNormFermion.size() - 1) outputNormalisation << ",";
  }
  outputNormalisation << "};" << std::endl;

  outputNormalisation << "}\n}\n}" << std::endl;

  outputNormalisation
      << "\n"
      << "#endif /* SRC_BARYO_CALCULATION_KFACTORS_GRID_KtildeInterpolation_H_ "
         "*/"
      << std::endl;

  outputNormalisation.close();
}

void CreateKfunctionsgrid()
{

  bool first = true;

  std::vector<std::vector<double>> K1p, K2p, K4p, K5p, K6p, K8p, K9p, K1m, K4m,
      K5m;

  std::vector<double> masslist, Tlist;
  for (double im = 0; im <= 5; im += 1e-3)
    masslist.push_back(im);
  for (double im = 10; im <= std::pow(200, 2); im += 5)
    masslist.push_back(im);
  for (double it = 10; it <= 250; it += 2)
    Tlist.push_back(it);

  for (auto im : masslist)
  {
    std::vector<double> K1p_row, K2p_row, K4p_row, K5p_row, K6p_row, K8p_row,
        K9p_row, K1m_row, K4m_row, K5m_row;
    for (auto it : Tlist)
    {
      K1p_row.push_back(Kfactors::K_integration(im, it, 1, 1));
      K2p_row.push_back(Kfactors::K_integration(im, it, 2, 1));
      K4p_row.push_back(Kfactors::K_integration(im, it, 4, 1));
      K5p_row.push_back(Kfactors::K_integration(im, it, 5, 1));
      K6p_row.push_back(Kfactors::K_integration(im, it, 6, 1));
      K8p_row.push_back(Kfactors::K_integration(im, it, 8, 1));
      K9p_row.push_back(Kfactors::K_integration(std::max(im, 1e-15), it, 9, 1));
      K1m_row.push_back(Kfactors::K_integration(im, it, 1, -1));
      K4m_row.push_back(Kfactors::K_integration(im, it, 4, -1));
      K5m_row.push_back(Kfactors::K_integration(im, it, 5, -1));
    }
    K1p.push_back(K1p_row);
    K2p.push_back(K2p_row);
    K4p.push_back(K4p_row);
    K5p.push_back(K5p_row);
    K6p.push_back(K6p_row);
    K8p.push_back(K8p_row);
    K9p.push_back(K9p_row);
    K1m.push_back(K1m_row);
    K4m.push_back(K4m_row);
    K5m.push_back(K5m_row);
  }

  std::ofstream header("Kfunctions_grid.h");
  header
      << "/*"
      << "\n"
      << " * Kfunctions_grid.h"
      << "\n"
      << " *"
      << "\n"
      << " *  Copyright (C) 2020  Philipp Basler, Margarete Mühlleitner and "
         "Jonas Müller"
      << "\n"
      << ""
      << "\n"
      << "      This program is free software: you can redistribute it and/or "
         "modify"
      << "\n"
      << "      it under the terms of the GNU General Public License as "
         "published by"
      << "\n"
      << "      the Free Software Foundation, either version 3 of the License, "
         "or"
      << "\n"
      << "      (at your option) any later version."
      << "\n"
      << ""
      << "\n"
      << "      This program is distributed in the hope that it will be useful,"
      << "\n"
      << "      but WITHOUT ANY WARRANTY; without even the implied warranty of"
      << "\n"
      << "      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the"
      << "\n"
      << "      GNU General Public License for more details."
      << "\n"
      << ""
      << "\n"
      << "      You should have received a copy of the GNU General Public "
         "License"
      << "\n"
      << "      along with this program.  If not, see "
         "<http://www.gnu.org/licenses/>."
      << "\n"
      << " */"
      << "\n";

  header << "#ifndef SRC_BARYO_CALCULATION_KFACTORS_GRID_KFUNCTIONS_GRID_H_NEW"
         << std::endl
         << "#define SRC_BARYO_CALCULATION_KFACTORS_GRID_KFUNCTIONS_GRID_H_NEW"
         << std::endl;
  header << "namespace BSMPT {" << std::endl
         << "namespace Kfactors {" << std::endl
         << "namespace Data {" << std::endl;

  header << "const int msg_size = " << masslist.size() << ";" << std::endl;
  header << "const int Tg_size = " << Tlist.size() << ";" << std::endl;
  header << "const double msg[" << masslist.size() << "] = {";
  first = true;
  for (const auto &x : masslist)
  {
    if (first)
    {
      first = false;
    }
    else
    {
      header << ",";
    }
    header << x;
  }
  header << "};" << std::endl;

  first = true;
  header << "const double Tg[" << Tlist.size() << "] = {";
  for (const auto &x : Tlist)
  {
    if (first)
    {
      first = false;
    }
    else
    {
      header << ",";
    }
    header << x;
  }
  header << "};" << std::endl;

  header << "const double K1p[" << masslist.size() << "][" << Tlist.size()
         << "] = {";
  for (std::size_t i = 0; i < masslist.size(); i++)
  {
    header << "{";
    for (std::size_t j = 0; j < Tlist.size(); j++)
    {
      header << K1p[i][j];
      if (j < Tlist.size() - 1) header << ",";
    }
    header << "}";
    if (i < masslist.size() - 1) header << ",";
  }
  header << "};" << std::endl;

  header << "const double K1m[" << masslist.size() << "][" << Tlist.size()
         << "] = {";
  for (std::size_t i = 0; i < masslist.size(); i++)
  {
    header << "{";
    for (std::size_t j = 0; j < Tlist.size(); j++)
    {
      header << K1m[i][j];
      if (j < Tlist.size() - 1) header << ",";
    }
    header << "}";
    if (i < masslist.size() - 1) header << ",";
  }
  header << "};" << std::endl;

  header << "const double K2p[" << masslist.size() << "][" << Tlist.size()
         << "] = {";
  for (std::size_t i = 0; i < masslist.size(); i++)
  {
    header << "{";
    for (std::size_t j = 0; j < Tlist.size(); j++)
    {
      header << K2p[i][j];
      if (j < Tlist.size() - 1) header << ",";
    }
    header << "}";
    if (i < masslist.size() - 1) header << ",";
  }
  header << "};" << std::endl;

  header << "const double K4p[" << masslist.size() << "][" << Tlist.size()
         << "] = {";
  for (std::size_t i = 0; i < masslist.size(); i++)
  {
    header << "{";
    for (std::size_t j = 0; j < Tlist.size(); j++)
    {
      header << K4p[i][j];
      if (j < Tlist.size() - 1) header << ",";
    }
    header << "}";
    if (i < masslist.size() - 1) header << ",";
  }
  header << "};" << std::endl;

  header << "const double K4m[" << masslist.size() << "][" << Tlist.size()
         << "] = {";
  for (std::size_t i = 0; i < masslist.size(); i++)
  {
    header << "{";
    for (std::size_t j = 0; j < Tlist.size(); j++)
    {
      header << K4m[i][j];
      if (j < Tlist.size() - 1) header << ",";
    }
    header << "}";
    if (i < masslist.size() - 1) header << ",";
  }
  header << "};" << std::endl;

  header << "const double K5p[" << masslist.size() << "][" << Tlist.size()
         << "] = {";
  for (std::size_t i = 0; i < masslist.size(); i++)
  {
    header << "{";
    for (std::size_t j = 0; j < Tlist.size(); j++)
    {
      header << K5p[i][j];
      if (j < Tlist.size() - 1) header << ",";
    }
    header << "}";
    if (i < masslist.size() - 1) header << ",";
  }
  header << "};" << std::endl;

  header << "const double K5m[" << masslist.size() << "][" << Tlist.size()
         << "] = {";
  for (std::size_t i = 0; i < masslist.size(); i++)
  {
    header << "{";
    for (std::size_t j = 0; j < Tlist.size(); j++)
    {
      header << K5m[i][j];
      if (j < Tlist.size() - 1) header << ",";
    }
    header << "}";
    if (i < masslist.size() - 1) header << ",";
  }
  header << "};" << std::endl;

  header << "const double K6p[" << masslist.size() << "][" << Tlist.size()
         << "] = {";
  for (std::size_t i = 0; i < masslist.size(); i++)
  {
    header << "{";
    for (std::size_t j = 0; j < Tlist.size(); j++)
    {
      header << K6p[i][j];
      if (j < Tlist.size() - 1) header << ",";
    }
    header << "}";
    if (i < masslist.size() - 1) header << ",";
  }
  header << "};" << std::endl;

  header << "const double K8p[" << masslist.size() << "][" << Tlist.size()
         << "] = {";
  for (std::size_t i = 0; i < masslist.size(); i++)
  {
    header << "{";
    for (std::size_t j = 0; j < Tlist.size(); j++)
    {
      header << K8p[i][j];
      if (j < Tlist.size() - 1) header << ",";
    }
    header << "}";
    if (i < masslist.size() - 1) header << ",";
  }
  header << "};" << std::endl;

  header << "const double K9p[" << masslist.size() << "][" << Tlist.size()
         << "] = {";
  for (std::size_t i = 0; i < masslist.size(); i++)
  {
    header << "{";
    for (std::size_t j = 0; j < Tlist.size(); j++)
    {
      header << K9p[i][j];
      if (j < Tlist.size() - 1) header << ",";
    }
    header << "}";
    if (i < masslist.size() - 1) header << ",";
  }
  header << "};" << std::endl;

  header << "}" << std::endl
         << "}" << std::endl
         << "}" << std::endl
         << "#endif" << std::endl;

  header.close();
}

int main()
try
{

  CreateKtildeInterpolationData();
  CreateKfunctionsgrid();

  return EXIT_SUCCESS;
}

catch (exception &e)
{
  std::cerr << e.what() << std::endl;
  return EXIT_FAILURE;
}
