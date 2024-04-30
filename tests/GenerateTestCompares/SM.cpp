// SPDX-FileCopyrightText: 2021 Philipp Basler
//
// SPDX-License-Identifier: GPL-3.0-or-later
#include "SM.h"
Compare_SM::Compare_SM()
{
  std::size_t NHiggs = 4;
  CheckTripleTree =
      Matrix3D{NHiggs, Matrix2D{NHiggs, std::vector<double>(NHiggs, 0)}};
  CheckTripleCW =
      Matrix3D{NHiggs, Matrix2D{NHiggs, std::vector<double>(NHiggs, 0)}};
  CheckTripleCT =
      Matrix3D{NHiggs, Matrix2D{NHiggs, std::vector<double>(NHiggs, 0)}};
  EWPTPerSetting[4].Tc = 159.082;
  EWPTPerSetting[4].vc = 22.6653;
  EWPTPerSetting[4].EWMinimum.push_back(-22.6653);
  EWPTPerSetting[1].Tc = 159.082;
  EWPTPerSetting[1].vc = 22.6654;
  EWPTPerSetting[1].EWMinimum.push_back(-22.6654);
  EWPTPerSetting[5].Tc = 159.082;
  EWPTPerSetting[5].vc = 22.6653;
  EWPTPerSetting[5].EWMinimum.push_back(-22.6653);
  EWPTPerSetting[2].Tc = 159.082;
  EWPTPerSetting[2].vc = 22.6654;
  EWPTPerSetting[2].EWMinimum.push_back(-22.6654);
  EWPTPerSetting[6].Tc = 159.082;
  EWPTPerSetting[6].vc = 22.6653;
  EWPTPerSetting[6].EWMinimum.push_back(-22.6653);
  EWPTPerSetting[3].Tc = 159.082;
  EWPTPerSetting[3].vc = 22.6654;
  EWPTPerSetting[3].EWMinimum.push_back(-22.6654);
  EWPTPerSetting[7].Tc = 159.082;
  EWPTPerSetting[7].vc = 22.6653;
  EWPTPerSetting[7].EWMinimum.push_back(-22.6653);
  CheckTripleTree.at(0).at(0).at(3) = 63.551;
  CheckTripleCT.at(0).at(0).at(3)   = -5.21625;
  CheckTripleCW.at(0).at(0).at(3)   = 5.21625;
  CheckTripleTree.at(0).at(3).at(0) = 63.551;
  CheckTripleCT.at(0).at(3).at(0)   = -5.21625;
  CheckTripleCW.at(0).at(3).at(0)   = 5.21625;
  CheckTripleTree.at(1).at(1).at(3) = 63.551;
  CheckTripleCT.at(1).at(1).at(3)   = -5.21625;
  CheckTripleCW.at(1).at(1).at(3)   = 5.21625;
  CheckTripleTree.at(1).at(3).at(1) = 63.551;
  CheckTripleCT.at(1).at(3).at(1)   = -5.21625;
  CheckTripleCW.at(1).at(3).at(1)   = 5.21625;
  CheckTripleTree.at(2).at(2).at(3) = 63.551;
  CheckTripleCT.at(2).at(2).at(3)   = -5.21625;
  CheckTripleCW.at(2).at(2).at(3)   = 5.21625;
  CheckTripleTree.at(2).at(3).at(2) = 63.551;
  CheckTripleCT.at(2).at(3).at(2)   = -5.21625;
  CheckTripleCW.at(2).at(3).at(2)   = 5.21625;
  CheckTripleTree.at(3).at(0).at(0) = 63.551;
  CheckTripleCT.at(3).at(0).at(0)   = -5.21625;
  CheckTripleCW.at(3).at(0).at(0)   = 5.21625;
  CheckTripleTree.at(3).at(1).at(1) = 63.551;
  CheckTripleCT.at(3).at(1).at(1)   = -5.21625;
  CheckTripleCW.at(3).at(1).at(1)   = 5.21625;
  CheckTripleTree.at(3).at(2).at(2) = 63.551;
  CheckTripleCT.at(3).at(2).at(2)   = -5.21625;
  CheckTripleCW.at(3).at(2).at(2)   = 5.21625;
  CheckTripleTree.at(3).at(3).at(3) = 190.653;
  CheckTripleCT.at(3).at(3).at(3)   = -15.6487;
  CheckTripleCW.at(3).at(3).at(3)   = -0.202437;
}
