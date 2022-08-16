// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
//
// SPDX-License-Identifier: GPL-3.0-or-later
#include <BSMPT/minimizer/Minimizer.h>
#include <map>
#include <vector>
class Compare_CPINTHEDARK
{
public:
  using Matrix3D = std::vector<std::vector<std::vector<double>>>;
  using Matrix2D = std::vector<std::vector<double>>;
  Compare_CPINTHEDARK();
  Matrix3D CheckTripleCT;
  Matrix3D CheckTripleCW;
  Matrix3D CheckTripleTree;
  std::map<int, BSMPT::Minimizer::EWPTReturnType> EWPTPerSetting;
};
