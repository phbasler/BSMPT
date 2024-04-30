// SPDX-FileCopyrightText: 2021 Philipp Basler
//
// SPDX-License-Identifier: GPL-3.0-or-later
#include <BSMPT/minimizer/Minimizer.h>
#include <map>
#include <vector>
class Compare_C2HDM
{
public:
  using Matrix3D = std::vector<std::vector<std::vector<double>>>;
  using Matrix2D = std::vector<std::vector<double>>;
  Compare_C2HDM();
  Matrix3D CheckTripleCT;
  Matrix3D CheckTripleCW;
  Matrix3D CheckTripleTree;
  std::map<int, BSMPT::Minimizer::EWPTReturnType> EWPTPerSetting;
  std::map<int, double> LWPerSetting;
  std::map<int, std::vector<double>> vevSymmetricPerSetting;
  std::map<int, std::vector<double>> etaPerSetting;
  const double testVW = 0.1;
};
