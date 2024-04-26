// SPDX-FileCopyrightText: 2021 Philipp Basler
//
// SPDX-License-Identifier: GPL-3.0-or-later
#include <BSMPT/minimizer/Minimizer.h>
#include <map>
#include <vector>
class Compare_CXSM
{
public:
  using Matrix3D = std::vector<std::vector<std::vector<double>>>;
  using Matrix2D = std::vector<std::vector<double>>;
  Compare_CXSM();
  Matrix3D CheckTripleCT;
  Matrix3D CheckTripleCW;
  Matrix3D CheckTripleTree;
  std::map<int, BSMPT::Minimizer::EWPTReturnType> EWPTPerSetting;
};
