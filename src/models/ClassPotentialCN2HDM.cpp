/*
 * ClassTemplate.cpp
 *
 *

		This program is free software: you can redistribute it and/or modify
		it under the terms of the GNU General Public License as published by
		the Free Software Foundation, either version 3 of the License, or
		(at your option) any later version.

		This program is distributed in the hope that it will be useful,
		but WITHOUT ANY WARRANTY; without even the implied warranty of
		MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
		GNU General Public License for more details.

		You should have received a copy of the GNU General Public License
		along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <BSMPT/models/ClassPotentialCN2HDM.h>
#include <BSMPT/models/IncludeAllModels.h>
using namespace Eigen;

/**
 * @file
 */

namespace BSMPT {
namespace Models {

/**
 * Here you have to adjust NNeutralHiggs, NChargedHiggs, nPar (number of Lagrangian parameters AFTER
 *  using the tadpole conditions),
 * nParCT (number of counterterms) as well as nVEV (number of VEVs for minimization)
 */
Class_Potential_CN2HDM::Class_Potential_CN2HDM ()
{
  Model = ModelID::ModelIDs::CN2HDM; // global int constant which will be used to tell the program which model is called
  NNeutralHiggs = 6; // number of neutral Higgs bosons at T = 0
  NChargedHiggs=4; // number of charged Higgs bosons  at T = 0 (all d.o.f.)

  nPar = 15;
  //nPar = 14; // number of parameters in the tree-Level Lagrangian


  nVEV=5; // number of VEVs to minimize the potential
  if(UseDMVEV) nVEV = 6;


  NHiggs = NNeutralHiggs+NChargedHiggs;

  nParCT = nPar+NHiggs; // number of parameters in the counterterm potential
  // nParCT = nPar+NHiggs+1;

  VevOrder.resize(nVEV);
  // Here you have to tell which scalar field gets which VEV.
  // ORDER: omega_1, omega_{CB}, omega_2, omega_{CP}, omega_s, {omega_{DM}}";
  VevOrder[0] = 2;
  VevOrder[1] = 4;
  VevOrder[2] = 6;
  VevOrder[3] = 7;
  VevOrder[4] = 8;
  if(UseDMVEV) VevOrder[5] = 9;

  // Set UseVTreeSimplified to use the tree-level potential defined in VTreeSimplified
  UseVTreeSimplified = false;

  // Set UseVCounterSimplified to use the counterterm potential defined in VCounterSimplified
  UseVCounterSimplified = false;

}

Class_Potential_CN2HDM::~Class_Potential_CN2HDM ()
{
}

/**
 * returns a string which tells the user the chronological order of the counterterms. Use this to
 * complement the legend of the given input file
 */
std::vector<std::string> Class_Potential_CN2HDM::addLegendCT() const{
    std::vector<std::string> labels;
    labels.push_back("dm112");
    labels.push_back("dm222");
    labels.push_back("dRem122");
    labels.push_back("dImm122");
    labels.push_back("dL1");
    labels.push_back("dL2");
    labels.push_back("dL3");
    labels.push_back("dL4");
    labels.push_back("dReL5");
    labels.push_back("dImL5");
    labels.push_back("dms2");
    labels.push_back("dL6");
    labels.push_back("dL7");
    labels.push_back("dL8");
    labels.push_back("dmDM2");
    labels.push_back("dT1");
    labels.push_back("dT2");
    labels.push_back("dT3");
    labels.push_back("dT4");
    labels.push_back("dT5");
    labels.push_back("dT6");
    labels.push_back("dT7");
    labels.push_back("dT8");
    labels.push_back("dT9");
    labels.push_back("dT10");
    return labels;
}

/**
 * returns a string which tells the user the chronological order of the VEVs and the critical temperature. Use this to
 * complement the legend of the given input file
 */
std::vector<std::string> Class_Potential_CN2HDM::addLegendTemp() const{
    std::vector<std::string> labels;
    labels.push_back("T_c");
    labels.push_back("v_c");
    labels.push_back("omega_c/T_c");
    labels.push_back("omega_1(T_c)");
    labels.push_back("omega_{CB}(T_c)");
    labels.push_back("omega_2(T_c)");
    labels.push_back("omega_{CP}(T_c)");
    labels.push_back("omega_s(T_c)");
    if(UseDMVEV) labels.push_back("omega_{DM}(T_c)");
    return  labels;
}

/**
 * returns a string which tells the user the chronological order of the Triple Higgs couplings. Use this to
 * complement the legend of the given input file
 *
 */
std::vector<std::string> Class_Potential_CN2HDM::addLegendTripleCouplings() const{
	std::vector<std::string> particles;
    std::vector<std::string> labels;
	particles.resize(NHiggs);
	//here you have to define the particle names in the vector particles


    particles[0] = "G^+";
    particles[1] = "G^-";
    particles[2] = "H^+";
    particles[3] = "H^-";
    particles[4] = "G^0";
    particles[5] = "h1";
    particles[6] = "h2";
    particles[7] = "h3";
    particles[8] = "h4";
    particles[9] = "hDM";

    for(size_t i=0;i<NHiggs;i++)
    {
        for(size_t j=i;j<NHiggs;j++)
        {
            for(size_t k=j;k<NHiggs;k++)
            {
                labels.push_back("Tree_"+particles.at(i)+particles.at(j)+particles.at(k));
                labels.push_back("CT_"+particles.at(i)+particles.at(j)+particles.at(k));
                labels.push_back("CW_"+particles.at(i)+particles.at(j)+particles.at(k));
            }
        }
    }

    return labels;
}

/**
 * returns a string which tells the user the chronological order of the VEVs. Use this to
 * complement the legend of the given input file
 */
std::vector<std::string> Class_Potential_CN2HDM::addLegendVEV() const{
    std::vector<std::string> labels;
	//out = "Your VEV order";
    labels.push_back("omega_1");
    labels.push_back("omega_{CB}");
    labels.push_back("omega_2");
    labels.push_back("omega_{CP}");
    labels.push_back("omega_s");
    if(UseDMVEV) labels.push_back("omega_{DM}");
    return labels;
}

/**
 * Reads the string linestr and sets the parameter point
 */
void Class_Potential_CN2HDM::ReadAndSet(const std::string& linestr, std::vector<double>& par )
{
	std::stringstream ss(linestr);
	double tmp;

    if (UseIndexCol){
        ss >> tmp;
    }
    // the column variables of the input file should be in this order. What comes in the columns after this point will be ignored. Order should be consistent and written in 'prepare_data_CN2HDM.py'

	for(int k=1;k<=14;k++)
    {
          ss>>tmp;
          if(k==1) par[0] = tmp; // Type
          else if(k==2) par[1] = tmp; // Re(m12^2)
          else if(k==3) par[2] = tmp; // L1
          else if(k==4) par[3] = tmp; // L2;
          else if(k==5) par[4] = tmp; // L3;
          else if(k==6) par[5] = tmp; // L4;
          else if(k==7) par[6] = tmp; // Re(L5)
          else if(k==8) par[7] = tmp; // Im(L5)
          else if(k==9) par[8] = tmp; // L6;
          else if(k==10) par[9] = tmp; // L7
          else if(k==11) par[10] = tmp; // L8
          else if(k==12) par[11] = tmp; // md^2
          else if(k==13) par[12] = tmp; // tan(beta)
          else if(k==14) par[13] = tmp; // v_s
    }


	set_gen(par); // This you have to call so that everything will be set
	return ;
}


/**
 * Set Class Object as well as the VEV configuration
 */
void Class_Potential_CN2HDM::set_gen(const std::vector<double>& par) {
    // define here also parameters that are not read in from the input parameters: angles; solutions of tadpole eqns etc.

    //independent parameters [0-13]: these that follow are defined in input file and read in.
    Type = static_cast<int>(par[0]);
    Rem122 = par[10];
    L1 = par[1];
    L2 = par[2];
    L3 = par[3];
    L4 = par[4];
    ReL5 = par[5];
    ImL5 = par[6];
    L6 = par[7];
    L7 = par[8];
    L8 = par[9];
    mds = par[11];
    TanBeta = par[12];
    vs = par[13];

    mDM2 = mds;

    // VEVs

    using std::sqrt;
    using std::pow;

    scale = C_vev0;

    vevTreeMin.resize(nVEV);
    vevTree.resize(NHiggs);

    C_CosBetaSquared = 1.0/(1+TanBeta*TanBeta);
    C_CosBeta = sqrt(C_CosBetaSquared);
    C_SinBetaSquared = (TanBeta*TanBeta)/(1+TanBeta*TanBeta);
    C_SinBeta = sqrt(C_SinBetaSquared);

    v1 = C_vev0 * C_CosBeta;
    v2 = C_vev0 * C_SinBeta;

    // Here you have to set the vector vevTreeMin. The vector vevTree will then be set by the function MinimizeOrderVEV

    // ORDER: omega_1, omega_{CB}, omega_2, omega_{CP}, omega_s, {omega_{DM}}";

    vevTreeMin[0] = v1;
    vevTreeMin[1] = 0;
    vevTreeMin[2] = v2;
    vevTreeMin[3] = 0;
    vevTreeMin[4] = vs;
    if(UseDMVEV) vevTreeMin[5] = 0;

    // Minimisation conditions

    Imm122 = v1 * v2 * ImL5 / 0.2e1;
    m112 = -0.1e1 / v1 * (0.2e1 * L1 * pow(v1, 0.3e1) + 0.2e1 * L3 * v1 * v2 * v2 + 0.2e1 * L4 * v2 * v2 * v1 + L7 * v1 * vs * vs + 0.2e1 * ReL5 * v1 * v2 * v2 - 0.4e1 * Rem122 * v2) / 0.4e1;
    m222 = -0.1e1 / v2 * (0.2e1 * L2 * pow(v2, 0.3e1) + 0.2e1 * L3 * v1 * v1 * v2 + 0.2e1 * L4 * v1 * v1 * v2 + L8 * v2 * vs * vs + 0.2e1 * ReL5 * v1 * v1 * v2 - 0.4e1 * Rem122 * v1) / 0.4e1;
    ms2 = -L6 * vs * vs / 0.4e1 - L7 * v1 * v1 / 0.2e1 - L8 * v2 * v2 / 0.2e1 + mDM2;


    vevTree=MinimizeOrderVEV(vevTreeMin);
	if(!SetCurvatureDone) SetCurvatureArrays();
}

/**
 * set your counterterm parameters from the entries of par as well as the entries of Curvature_Higgs_CT_L1 to
 * Curvature_Higgs_CT_L4.
 */
void Class_Potential_CN2HDM::set_CT_Pot_Par(const std::vector<double>& par){

    dm112 = par[0];
    dm222 = par[1];
    dRem122 = par[2];
    dImm122 = par[3];
    dL1 = par[4];
    dL2 = par[5];
    dL3 = par[6];
    dL4 = par[7];
    dReL5 = par[8];
    dImL5 = par[9];
    dms2 = par[10];
    dL6 = par[11];
    dL7 = par[12];
    dL8 = par[13];
    dmDM2 = par[14];
    dT1 = par[15];
    dT2 = par[16];
    dT3 = par[17];
    dT4 = par[18];
    dT5 = par[19];
    dT6 = par[20];
    dT7 = par[21];
    dT8 = par[22];
    dT9 = par[23];
    dT10 = par[24];

    Curvature_Higgs_CT_L1[0] = dT1;
    Curvature_Higgs_CT_L1[1] = dT2;
    Curvature_Higgs_CT_L1[2] = dT3;
    Curvature_Higgs_CT_L1[3] = dT4;
    Curvature_Higgs_CT_L1[4] = dT5;
    Curvature_Higgs_CT_L1[5] = dT6;
    Curvature_Higgs_CT_L1[6] = dT7;
    Curvature_Higgs_CT_L1[7] = dT8;
    Curvature_Higgs_CT_L1[8] = dT9;
    Curvature_Higgs_CT_L1[9] = dT10;

    Curvature_Higgs_CT_L2[0][0] = dm112;
    Curvature_Higgs_CT_L2[0][5] = dImm122;
    Curvature_Higgs_CT_L2[1][1] = dm112;
    Curvature_Higgs_CT_L2[1][4] = -dImm122;
    Curvature_Higgs_CT_L2[1][5] = -dRem122;
    Curvature_Higgs_CT_L2[2][2] = dm112;
    Curvature_Higgs_CT_L2[2][6] = -dRem122;
    Curvature_Higgs_CT_L2[2][7] = dImm122;
    Curvature_Higgs_CT_L2[3][3] = dm112;
    Curvature_Higgs_CT_L2[3][6] = -dImm122;
    Curvature_Higgs_CT_L2[3][7] = -dRem122;
    Curvature_Higgs_CT_L2[4][0] = -dRem122;
    Curvature_Higgs_CT_L2[4][1] = -dImm122;
    Curvature_Higgs_CT_L2[4][4] = dm222;
    Curvature_Higgs_CT_L2[5][0] = dImm122;
    Curvature_Higgs_CT_L2[5][1] = -dRem122;
    Curvature_Higgs_CT_L2[5][5] = dm222;
    Curvature_Higgs_CT_L2[6][2] = -dRem122;
    Curvature_Higgs_CT_L2[6][3] = -dImm122;
    Curvature_Higgs_CT_L2[6][6] = dm222;
    Curvature_Higgs_CT_L2[7][2] = dImm122;
    Curvature_Higgs_CT_L2[7][3] = -dRem122;
    Curvature_Higgs_CT_L2[7][7] = dm222;
    Curvature_Higgs_CT_L2[8][8] = dms2 / 0.2e1 - dmDM2 / 0.2e1;
    Curvature_Higgs_CT_L2[9][9] = dms2 / 0.2e1 + dmDM2 / 0.2e1;

    Curvature_Higgs_CT_L4[0][0][0][0] = 3 * dL1;
    Curvature_Higgs_CT_L4[0][0][1][1] = dL1;
    Curvature_Higgs_CT_L4[0][0][2][2] = dL1;
    Curvature_Higgs_CT_L4[0][0][3][3] = dL1;
    Curvature_Higgs_CT_L4[0][0][4][4] = dL3 + dL4 + dReL5;
    Curvature_Higgs_CT_L4[0][0][4][5] = -dImL5;
    Curvature_Higgs_CT_L4[0][0][5][4] = -dImL5;
    Curvature_Higgs_CT_L4[0][0][5][5] = dL3 + dL4 - dReL5;
    Curvature_Higgs_CT_L4[0][0][6][6] = dL3;
    Curvature_Higgs_CT_L4[0][0][7][7] = dL3;
    Curvature_Higgs_CT_L4[0][0][8][8] = dL7 / 0.2e1;
    Curvature_Higgs_CT_L4[0][0][9][9] = dL7 / 0.2e1;
    Curvature_Higgs_CT_L4[0][1][0][1] = dL1;
    Curvature_Higgs_CT_L4[0][1][1][0] = dL1;
    Curvature_Higgs_CT_L4[0][1][4][4] = dImL5;
    Curvature_Higgs_CT_L4[0][1][4][5] = dReL5;
    Curvature_Higgs_CT_L4[0][1][5][4] = dReL5;
    Curvature_Higgs_CT_L4[0][1][5][5] = -dImL5;
    Curvature_Higgs_CT_L4[0][2][0][2] = dL1;
    Curvature_Higgs_CT_L4[0][2][2][0] = dL1;
    Curvature_Higgs_CT_L4[0][2][4][6] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[0][2][4][7] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[0][2][5][6] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[0][2][5][7] = dL4 / 0.2e1 - dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[0][2][6][4] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[0][2][6][5] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[0][2][7][4] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[0][2][7][5] = dL4 / 0.2e1 - dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[0][3][0][3] = dL1;
    Curvature_Higgs_CT_L4[0][3][3][0] = dL1;
    Curvature_Higgs_CT_L4[0][3][4][6] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[0][3][4][7] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[0][3][5][6] = -dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[0][3][5][7] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[0][3][6][4] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[0][3][6][5] = -dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[0][3][7][4] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[0][3][7][5] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[0][4][0][4] = dL3 + dL4 + dReL5;
    Curvature_Higgs_CT_L4[0][4][0][5] = -dImL5;
    Curvature_Higgs_CT_L4[0][4][1][4] = dImL5;
    Curvature_Higgs_CT_L4[0][4][1][5] = dReL5;
    Curvature_Higgs_CT_L4[0][4][2][6] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[0][4][2][7] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[0][4][3][6] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[0][4][3][7] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[0][4][4][0] = dL3 + dL4 + dReL5;
    Curvature_Higgs_CT_L4[0][4][4][1] = dImL5;
    Curvature_Higgs_CT_L4[0][4][5][0] = -dImL5;
    Curvature_Higgs_CT_L4[0][4][5][1] = dReL5;
    Curvature_Higgs_CT_L4[0][4][6][2] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[0][4][6][3] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[0][4][7][2] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[0][4][7][3] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[0][5][0][4] = -dImL5;
    Curvature_Higgs_CT_L4[0][5][0][5] = dL3 + dL4 - dReL5;
    Curvature_Higgs_CT_L4[0][5][1][4] = dReL5;
    Curvature_Higgs_CT_L4[0][5][1][5] = -dImL5;
    Curvature_Higgs_CT_L4[0][5][2][6] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[0][5][2][7] = dL4 / 0.2e1 - dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[0][5][3][6] = -dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[0][5][3][7] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[0][5][4][0] = -dImL5;
    Curvature_Higgs_CT_L4[0][5][4][1] = dReL5;
    Curvature_Higgs_CT_L4[0][5][5][0] = dL3 + dL4 - dReL5;
    Curvature_Higgs_CT_L4[0][5][5][1] = -dImL5;
    Curvature_Higgs_CT_L4[0][5][6][2] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[0][5][6][3] = -dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[0][5][7][2] = dL4 / 0.2e1 - dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[0][5][7][3] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[0][6][0][6] = dL3;
    Curvature_Higgs_CT_L4[0][6][2][4] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[0][6][2][5] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[0][6][3][4] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[0][6][3][5] = -dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[0][6][4][2] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[0][6][4][3] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[0][6][5][2] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[0][6][5][3] = -dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[0][6][6][0] = dL3;
    Curvature_Higgs_CT_L4[0][7][0][7] = dL3;
    Curvature_Higgs_CT_L4[0][7][2][4] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[0][7][2][5] = dL4 / 0.2e1 - dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[0][7][3][4] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[0][7][3][5] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[0][7][4][2] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[0][7][4][3] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[0][7][5][2] = dL4 / 0.2e1 - dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[0][7][5][3] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[0][7][7][0] = dL3;
    Curvature_Higgs_CT_L4[0][8][0][8] = dL7 / 0.2e1;
    Curvature_Higgs_CT_L4[0][8][8][0] = dL7 / 0.2e1;
    Curvature_Higgs_CT_L4[0][9][0][9] = dL7 / 0.2e1;
    Curvature_Higgs_CT_L4[0][9][9][0] = dL7 / 0.2e1;
    Curvature_Higgs_CT_L4[1][0][0][1] = dL1;
    Curvature_Higgs_CT_L4[1][0][1][0] = dL1;
    Curvature_Higgs_CT_L4[1][0][4][4] = dImL5;
    Curvature_Higgs_CT_L4[1][0][4][5] = dReL5;
    Curvature_Higgs_CT_L4[1][0][5][4] = dReL5;
    Curvature_Higgs_CT_L4[1][0][5][5] = -dImL5;
    Curvature_Higgs_CT_L4[1][1][0][0] = dL1;
    Curvature_Higgs_CT_L4[1][1][1][1] = 3 * dL1;
    Curvature_Higgs_CT_L4[1][1][2][2] = dL1;
    Curvature_Higgs_CT_L4[1][1][3][3] = dL1;
    Curvature_Higgs_CT_L4[1][1][4][4] = dL3 + dL4 - dReL5;
    Curvature_Higgs_CT_L4[1][1][4][5] = dImL5;
    Curvature_Higgs_CT_L4[1][1][5][4] = dImL5;
    Curvature_Higgs_CT_L4[1][1][5][5] = dL3 + dL4 + dReL5;
    Curvature_Higgs_CT_L4[1][1][6][6] = dL3;
    Curvature_Higgs_CT_L4[1][1][7][7] = dL3;
    Curvature_Higgs_CT_L4[1][1][8][8] = dL7 / 0.2e1;
    Curvature_Higgs_CT_L4[1][1][9][9] = dL7 / 0.2e1;
    Curvature_Higgs_CT_L4[1][2][1][2] = dL1;
    Curvature_Higgs_CT_L4[1][2][2][1] = dL1;
    Curvature_Higgs_CT_L4[1][2][4][6] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[1][2][4][7] = -dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[1][2][5][6] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[1][2][5][7] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[1][2][6][4] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[1][2][6][5] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[1][2][7][4] = -dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[1][2][7][5] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[1][3][1][3] = dL1;
    Curvature_Higgs_CT_L4[1][3][3][1] = dL1;
    Curvature_Higgs_CT_L4[1][3][4][6] = dL4 / 0.2e1 - dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[1][3][4][7] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[1][3][5][6] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[1][3][5][7] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[1][3][6][4] = dL4 / 0.2e1 - dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[1][3][6][5] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[1][3][7][4] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[1][3][7][5] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[1][4][0][4] = dImL5;
    Curvature_Higgs_CT_L4[1][4][0][5] = dReL5;
    Curvature_Higgs_CT_L4[1][4][1][4] = dL3 + dL4 - dReL5;
    Curvature_Higgs_CT_L4[1][4][1][5] = dImL5;
    Curvature_Higgs_CT_L4[1][4][2][6] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[1][4][2][7] = -dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[1][4][3][6] = dL4 / 0.2e1 - dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[1][4][3][7] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[1][4][4][0] = dImL5;
    Curvature_Higgs_CT_L4[1][4][4][1] = dL3 + dL4 - dReL5;
    Curvature_Higgs_CT_L4[1][4][5][0] = dReL5;
    Curvature_Higgs_CT_L4[1][4][5][1] = dImL5;
    Curvature_Higgs_CT_L4[1][4][6][2] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[1][4][6][3] = dL4 / 0.2e1 - dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[1][4][7][2] = -dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[1][4][7][3] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[1][5][0][4] = dReL5;
    Curvature_Higgs_CT_L4[1][5][0][5] = -dImL5;
    Curvature_Higgs_CT_L4[1][5][1][4] = dImL5;
    Curvature_Higgs_CT_L4[1][5][1][5] = dL3 + dL4 + dReL5;
    Curvature_Higgs_CT_L4[1][5][2][6] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[1][5][2][7] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[1][5][3][6] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[1][5][3][7] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[1][5][4][0] = dReL5;
    Curvature_Higgs_CT_L4[1][5][4][1] = dImL5;
    Curvature_Higgs_CT_L4[1][5][5][0] = -dImL5;
    Curvature_Higgs_CT_L4[1][5][5][1] = dL3 + dL4 + dReL5;
    Curvature_Higgs_CT_L4[1][5][6][2] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[1][5][6][3] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[1][5][7][2] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[1][5][7][3] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[1][6][1][6] = dL3;
    Curvature_Higgs_CT_L4[1][6][2][4] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[1][6][2][5] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[1][6][3][4] = dL4 / 0.2e1 - dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[1][6][3][5] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[1][6][4][2] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[1][6][4][3] = dL4 / 0.2e1 - dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[1][6][5][2] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[1][6][5][3] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[1][6][6][1] = dL3;
    Curvature_Higgs_CT_L4[1][7][1][7] = dL3;
    Curvature_Higgs_CT_L4[1][7][2][4] = -dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[1][7][2][5] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[1][7][3][4] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[1][7][3][5] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[1][7][4][2] = -dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[1][7][4][3] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[1][7][5][2] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[1][7][5][3] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[1][7][7][1] = dL3;
    Curvature_Higgs_CT_L4[1][8][1][8] = dL7 / 0.2e1;
    Curvature_Higgs_CT_L4[1][8][8][1] = dL7 / 0.2e1;
    Curvature_Higgs_CT_L4[1][9][1][9] = dL7 / 0.2e1;
    Curvature_Higgs_CT_L4[1][9][9][1] = dL7 / 0.2e1;
    Curvature_Higgs_CT_L4[2][0][0][2] = dL1;
    Curvature_Higgs_CT_L4[2][0][2][0] = dL1;
    Curvature_Higgs_CT_L4[2][0][4][6] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[2][0][4][7] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[2][0][5][6] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[2][0][5][7] = dL4 / 0.2e1 - dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[2][0][6][4] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[2][0][6][5] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[2][0][7][4] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[2][0][7][5] = dL4 / 0.2e1 - dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[2][1][1][2] = dL1;
    Curvature_Higgs_CT_L4[2][1][2][1] = dL1;
    Curvature_Higgs_CT_L4[2][1][4][6] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[2][1][4][7] = -dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[2][1][5][6] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[2][1][5][7] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[2][1][6][4] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[2][1][6][5] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[2][1][7][4] = -dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[2][1][7][5] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[2][2][0][0] = dL1;
    Curvature_Higgs_CT_L4[2][2][1][1] = dL1;
    Curvature_Higgs_CT_L4[2][2][2][2] = 3 * dL1;
    Curvature_Higgs_CT_L4[2][2][3][3] = dL1;
    Curvature_Higgs_CT_L4[2][2][4][4] = dL3;
    Curvature_Higgs_CT_L4[2][2][5][5] = dL3;
    Curvature_Higgs_CT_L4[2][2][6][6] = dL3 + dL4 + dReL5;
    Curvature_Higgs_CT_L4[2][2][6][7] = -dImL5;
    Curvature_Higgs_CT_L4[2][2][7][6] = -dImL5;
    Curvature_Higgs_CT_L4[2][2][7][7] = dL3 + dL4 - dReL5;
    Curvature_Higgs_CT_L4[2][2][8][8] = dL7 / 0.2e1;
    Curvature_Higgs_CT_L4[2][2][9][9] = dL7 / 0.2e1;
    Curvature_Higgs_CT_L4[2][3][2][3] = dL1;
    Curvature_Higgs_CT_L4[2][3][3][2] = dL1;
    Curvature_Higgs_CT_L4[2][3][6][6] = dImL5;
    Curvature_Higgs_CT_L4[2][3][6][7] = dReL5;
    Curvature_Higgs_CT_L4[2][3][7][6] = dReL5;
    Curvature_Higgs_CT_L4[2][3][7][7] = -dImL5;
    Curvature_Higgs_CT_L4[2][4][0][6] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[2][4][0][7] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[2][4][1][6] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[2][4][1][7] = -dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[2][4][2][4] = dL3;
    Curvature_Higgs_CT_L4[2][4][4][2] = dL3;
    Curvature_Higgs_CT_L4[2][4][6][0] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[2][4][6][1] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[2][4][7][0] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[2][4][7][1] = -dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[2][5][0][6] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[2][5][0][7] = dL4 / 0.2e1 - dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[2][5][1][6] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[2][5][1][7] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[2][5][2][5] = dL3;
    Curvature_Higgs_CT_L4[2][5][5][2] = dL3;
    Curvature_Higgs_CT_L4[2][5][6][0] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[2][5][6][1] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[2][5][7][0] = dL4 / 0.2e1 - dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[2][5][7][1] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[2][6][0][4] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[2][6][0][5] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[2][6][1][4] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[2][6][1][5] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[2][6][2][6] = dL3 + dL4 + dReL5;
    Curvature_Higgs_CT_L4[2][6][2][7] = -dImL5;
    Curvature_Higgs_CT_L4[2][6][3][6] = dImL5;
    Curvature_Higgs_CT_L4[2][6][3][7] = dReL5;
    Curvature_Higgs_CT_L4[2][6][4][0] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[2][6][4][1] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[2][6][5][0] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[2][6][5][1] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[2][6][6][2] = dL3 + dL4 + dReL5;
    Curvature_Higgs_CT_L4[2][6][6][3] = dImL5;
    Curvature_Higgs_CT_L4[2][6][7][2] = -dImL5;
    Curvature_Higgs_CT_L4[2][6][7][3] = dReL5;
    Curvature_Higgs_CT_L4[2][7][0][4] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[2][7][0][5] = dL4 / 0.2e1 - dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[2][7][1][4] = -dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[2][7][1][5] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[2][7][2][6] = -dImL5;
    Curvature_Higgs_CT_L4[2][7][2][7] = dL3 + dL4 - dReL5;
    Curvature_Higgs_CT_L4[2][7][3][6] = dReL5;
    Curvature_Higgs_CT_L4[2][7][3][7] = -dImL5;
    Curvature_Higgs_CT_L4[2][7][4][0] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[2][7][4][1] = -dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[2][7][5][0] = dL4 / 0.2e1 - dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[2][7][5][1] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[2][7][6][2] = -dImL5;
    Curvature_Higgs_CT_L4[2][7][6][3] = dReL5;
    Curvature_Higgs_CT_L4[2][7][7][2] = dL3 + dL4 - dReL5;
    Curvature_Higgs_CT_L4[2][7][7][3] = -dImL5;
    Curvature_Higgs_CT_L4[2][8][2][8] = dL7 / 0.2e1;
    Curvature_Higgs_CT_L4[2][8][8][2] = dL7 / 0.2e1;
    Curvature_Higgs_CT_L4[2][9][2][9] = dL7 / 0.2e1;
    Curvature_Higgs_CT_L4[2][9][9][2] = dL7 / 0.2e1;
    Curvature_Higgs_CT_L4[3][0][0][3] = dL1;
    Curvature_Higgs_CT_L4[3][0][3][0] = dL1;
    Curvature_Higgs_CT_L4[3][0][4][6] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[3][0][4][7] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[3][0][5][6] = -dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[3][0][5][7] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[3][0][6][4] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[3][0][6][5] = -dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[3][0][7][4] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[3][0][7][5] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[3][1][1][3] = dL1;
    Curvature_Higgs_CT_L4[3][1][3][1] = dL1;
    Curvature_Higgs_CT_L4[3][1][4][6] = dL4 / 0.2e1 - dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[3][1][4][7] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[3][1][5][6] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[3][1][5][7] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[3][1][6][4] = dL4 / 0.2e1 - dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[3][1][6][5] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[3][1][7][4] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[3][1][7][5] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[3][2][2][3] = dL1;
    Curvature_Higgs_CT_L4[3][2][3][2] = dL1;
    Curvature_Higgs_CT_L4[3][2][6][6] = dImL5;
    Curvature_Higgs_CT_L4[3][2][6][7] = dReL5;
    Curvature_Higgs_CT_L4[3][2][7][6] = dReL5;
    Curvature_Higgs_CT_L4[3][2][7][7] = -dImL5;
    Curvature_Higgs_CT_L4[3][3][0][0] = dL1;
    Curvature_Higgs_CT_L4[3][3][1][1] = dL1;
    Curvature_Higgs_CT_L4[3][3][2][2] = dL1;
    Curvature_Higgs_CT_L4[3][3][3][3] = 3 * dL1;
    Curvature_Higgs_CT_L4[3][3][4][4] = dL3;
    Curvature_Higgs_CT_L4[3][3][5][5] = dL3;
    Curvature_Higgs_CT_L4[3][3][6][6] = dL3 + dL4 - dReL5;
    Curvature_Higgs_CT_L4[3][3][6][7] = dImL5;
    Curvature_Higgs_CT_L4[3][3][7][6] = dImL5;
    Curvature_Higgs_CT_L4[3][3][7][7] = dL3 + dL4 + dReL5;
    Curvature_Higgs_CT_L4[3][3][8][8] = dL7 / 0.2e1;
    Curvature_Higgs_CT_L4[3][3][9][9] = dL7 / 0.2e1;
    Curvature_Higgs_CT_L4[3][4][0][6] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[3][4][0][7] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[3][4][1][6] = dL4 / 0.2e1 - dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[3][4][1][7] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[3][4][3][4] = dL3;
    Curvature_Higgs_CT_L4[3][4][4][3] = dL3;
    Curvature_Higgs_CT_L4[3][4][6][0] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[3][4][6][1] = dL4 / 0.2e1 - dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[3][4][7][0] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[3][4][7][1] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[3][5][0][6] = -dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[3][5][0][7] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[3][5][1][6] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[3][5][1][7] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[3][5][3][5] = dL3;
    Curvature_Higgs_CT_L4[3][5][5][3] = dL3;
    Curvature_Higgs_CT_L4[3][5][6][0] = -dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[3][5][6][1] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[3][5][7][0] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[3][5][7][1] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[3][6][0][4] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[3][6][0][5] = -dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[3][6][1][4] = dL4 / 0.2e1 - dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[3][6][1][5] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[3][6][2][6] = dImL5;
    Curvature_Higgs_CT_L4[3][6][2][7] = dReL5;
    Curvature_Higgs_CT_L4[3][6][3][6] = dL3 + dL4 - dReL5;
    Curvature_Higgs_CT_L4[3][6][3][7] = dImL5;
    Curvature_Higgs_CT_L4[3][6][4][0] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[3][6][4][1] = dL4 / 0.2e1 - dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[3][6][5][0] = -dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[3][6][5][1] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[3][6][6][2] = dImL5;
    Curvature_Higgs_CT_L4[3][6][6][3] = dL3 + dL4 - dReL5;
    Curvature_Higgs_CT_L4[3][6][7][2] = dReL5;
    Curvature_Higgs_CT_L4[3][6][7][3] = dImL5;
    Curvature_Higgs_CT_L4[3][7][0][4] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[3][7][0][5] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[3][7][1][4] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[3][7][1][5] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[3][7][2][6] = dReL5;
    Curvature_Higgs_CT_L4[3][7][2][7] = -dImL5;
    Curvature_Higgs_CT_L4[3][7][3][6] = dImL5;
    Curvature_Higgs_CT_L4[3][7][3][7] = dL3 + dL4 + dReL5;
    Curvature_Higgs_CT_L4[3][7][4][0] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[3][7][4][1] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[3][7][5][0] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[3][7][5][1] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[3][7][6][2] = dReL5;
    Curvature_Higgs_CT_L4[3][7][6][3] = dImL5;
    Curvature_Higgs_CT_L4[3][7][7][2] = -dImL5;
    Curvature_Higgs_CT_L4[3][7][7][3] = dL3 + dL4 + dReL5;
    Curvature_Higgs_CT_L4[3][8][3][8] = dL7 / 0.2e1;
    Curvature_Higgs_CT_L4[3][8][8][3] = dL7 / 0.2e1;
    Curvature_Higgs_CT_L4[3][9][3][9] = dL7 / 0.2e1;
    Curvature_Higgs_CT_L4[3][9][9][3] = dL7 / 0.2e1;
    Curvature_Higgs_CT_L4[4][0][0][4] = dL3 + dL4 + dReL5;
    Curvature_Higgs_CT_L4[4][0][0][5] = -dImL5;
    Curvature_Higgs_CT_L4[4][0][1][4] = dImL5;
    Curvature_Higgs_CT_L4[4][0][1][5] = dReL5;
    Curvature_Higgs_CT_L4[4][0][2][6] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[4][0][2][7] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[4][0][3][6] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[4][0][3][7] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[4][0][4][0] = dL3 + dL4 + dReL5;
    Curvature_Higgs_CT_L4[4][0][4][1] = dImL5;
    Curvature_Higgs_CT_L4[4][0][5][0] = -dImL5;
    Curvature_Higgs_CT_L4[4][0][5][1] = dReL5;
    Curvature_Higgs_CT_L4[4][0][6][2] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[4][0][6][3] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[4][0][7][2] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[4][0][7][3] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[4][1][0][4] = dImL5;
    Curvature_Higgs_CT_L4[4][1][0][5] = dReL5;
    Curvature_Higgs_CT_L4[4][1][1][4] = dL3 + dL4 - dReL5;
    Curvature_Higgs_CT_L4[4][1][1][5] = dImL5;
    Curvature_Higgs_CT_L4[4][1][2][6] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[4][1][2][7] = -dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[4][1][3][6] = dL4 / 0.2e1 - dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[4][1][3][7] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[4][1][4][0] = dImL5;
    Curvature_Higgs_CT_L4[4][1][4][1] = dL3 + dL4 - dReL5;
    Curvature_Higgs_CT_L4[4][1][5][0] = dReL5;
    Curvature_Higgs_CT_L4[4][1][5][1] = dImL5;
    Curvature_Higgs_CT_L4[4][1][6][2] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[4][1][6][3] = dL4 / 0.2e1 - dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[4][1][7][2] = -dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[4][1][7][3] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[4][2][0][6] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[4][2][0][7] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[4][2][1][6] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[4][2][1][7] = -dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[4][2][2][4] = dL3;
    Curvature_Higgs_CT_L4[4][2][4][2] = dL3;
    Curvature_Higgs_CT_L4[4][2][6][0] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[4][2][6][1] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[4][2][7][0] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[4][2][7][1] = -dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[4][3][0][6] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[4][3][0][7] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[4][3][1][6] = dL4 / 0.2e1 - dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[4][3][1][7] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[4][3][3][4] = dL3;
    Curvature_Higgs_CT_L4[4][3][4][3] = dL3;
    Curvature_Higgs_CT_L4[4][3][6][0] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[4][3][6][1] = dL4 / 0.2e1 - dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[4][3][7][0] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[4][3][7][1] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[4][4][0][0] = dL3 + dL4 + dReL5;
    Curvature_Higgs_CT_L4[4][4][0][1] = dImL5;
    Curvature_Higgs_CT_L4[4][4][1][0] = dImL5;
    Curvature_Higgs_CT_L4[4][4][1][1] = dL3 + dL4 - dReL5;
    Curvature_Higgs_CT_L4[4][4][2][2] = dL3;
    Curvature_Higgs_CT_L4[4][4][3][3] = dL3;
    Curvature_Higgs_CT_L4[4][4][4][4] = 3 * dL2;
    Curvature_Higgs_CT_L4[4][4][5][5] = dL2;
    Curvature_Higgs_CT_L4[4][4][6][6] = dL2;
    Curvature_Higgs_CT_L4[4][4][7][7] = dL2;
    Curvature_Higgs_CT_L4[4][4][8][8] = dL8 / 0.2e1;
    Curvature_Higgs_CT_L4[4][4][9][9] = dL8 / 0.2e1;
    Curvature_Higgs_CT_L4[4][5][0][0] = -dImL5;
    Curvature_Higgs_CT_L4[4][5][0][1] = dReL5;
    Curvature_Higgs_CT_L4[4][5][1][0] = dReL5;
    Curvature_Higgs_CT_L4[4][5][1][1] = dImL5;
    Curvature_Higgs_CT_L4[4][5][4][5] = dL2;
    Curvature_Higgs_CT_L4[4][5][5][4] = dL2;
    Curvature_Higgs_CT_L4[4][6][0][2] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[4][6][0][3] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[4][6][1][2] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[4][6][1][3] = dL4 / 0.2e1 - dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[4][6][2][0] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[4][6][2][1] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[4][6][3][0] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[4][6][3][1] = dL4 / 0.2e1 - dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[4][6][4][6] = dL2;
    Curvature_Higgs_CT_L4[4][6][6][4] = dL2;
    Curvature_Higgs_CT_L4[4][7][0][2] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[4][7][0][3] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[4][7][1][2] = -dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[4][7][1][3] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[4][7][2][0] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[4][7][2][1] = -dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[4][7][3][0] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[4][7][3][1] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[4][7][4][7] = dL2;
    Curvature_Higgs_CT_L4[4][7][7][4] = dL2;
    Curvature_Higgs_CT_L4[4][8][4][8] = dL8 / 0.2e1;
    Curvature_Higgs_CT_L4[4][8][8][4] = dL8 / 0.2e1;
    Curvature_Higgs_CT_L4[4][9][4][9] = dL8 / 0.2e1;
    Curvature_Higgs_CT_L4[4][9][9][4] = dL8 / 0.2e1;
    Curvature_Higgs_CT_L4[5][0][0][4] = -dImL5;
    Curvature_Higgs_CT_L4[5][0][0][5] = dL3 + dL4 - dReL5;
    Curvature_Higgs_CT_L4[5][0][1][4] = dReL5;
    Curvature_Higgs_CT_L4[5][0][1][5] = -dImL5;
    Curvature_Higgs_CT_L4[5][0][2][6] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[5][0][2][7] = dL4 / 0.2e1 - dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[5][0][3][6] = -dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[5][0][3][7] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[5][0][4][0] = -dImL5;
    Curvature_Higgs_CT_L4[5][0][4][1] = dReL5;
    Curvature_Higgs_CT_L4[5][0][5][0] = dL3 + dL4 - dReL5;
    Curvature_Higgs_CT_L4[5][0][5][1] = -dImL5;
    Curvature_Higgs_CT_L4[5][0][6][2] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[5][0][6][3] = -dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[5][0][7][2] = dL4 / 0.2e1 - dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[5][0][7][3] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[5][1][0][4] = dReL5;
    Curvature_Higgs_CT_L4[5][1][0][5] = -dImL5;
    Curvature_Higgs_CT_L4[5][1][1][4] = dImL5;
    Curvature_Higgs_CT_L4[5][1][1][5] = dL3 + dL4 + dReL5;
    Curvature_Higgs_CT_L4[5][1][2][6] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[5][1][2][7] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[5][1][3][6] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[5][1][3][7] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[5][1][4][0] = dReL5;
    Curvature_Higgs_CT_L4[5][1][4][1] = dImL5;
    Curvature_Higgs_CT_L4[5][1][5][0] = -dImL5;
    Curvature_Higgs_CT_L4[5][1][5][1] = dL3 + dL4 + dReL5;
    Curvature_Higgs_CT_L4[5][1][6][2] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[5][1][6][3] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[5][1][7][2] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[5][1][7][3] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[5][2][0][6] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[5][2][0][7] = dL4 / 0.2e1 - dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[5][2][1][6] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[5][2][1][7] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[5][2][2][5] = dL3;
    Curvature_Higgs_CT_L4[5][2][5][2] = dL3;
    Curvature_Higgs_CT_L4[5][2][6][0] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[5][2][6][1] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[5][2][7][0] = dL4 / 0.2e1 - dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[5][2][7][1] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[5][3][0][6] = -dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[5][3][0][7] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[5][3][1][6] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[5][3][1][7] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[5][3][3][5] = dL3;
    Curvature_Higgs_CT_L4[5][3][5][3] = dL3;
    Curvature_Higgs_CT_L4[5][3][6][0] = -dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[5][3][6][1] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[5][3][7][0] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[5][3][7][1] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[5][4][0][0] = -dImL5;
    Curvature_Higgs_CT_L4[5][4][0][1] = dReL5;
    Curvature_Higgs_CT_L4[5][4][1][0] = dReL5;
    Curvature_Higgs_CT_L4[5][4][1][1] = dImL5;
    Curvature_Higgs_CT_L4[5][4][4][5] = dL2;
    Curvature_Higgs_CT_L4[5][4][5][4] = dL2;
    Curvature_Higgs_CT_L4[5][5][0][0] = dL3 + dL4 - dReL5;
    Curvature_Higgs_CT_L4[5][5][0][1] = -dImL5;
    Curvature_Higgs_CT_L4[5][5][1][0] = -dImL5;
    Curvature_Higgs_CT_L4[5][5][1][1] = dL3 + dL4 + dReL5;
    Curvature_Higgs_CT_L4[5][5][2][2] = dL3;
    Curvature_Higgs_CT_L4[5][5][3][3] = dL3;
    Curvature_Higgs_CT_L4[5][5][4][4] = dL2;
    Curvature_Higgs_CT_L4[5][5][5][5] = 3 * dL2;
    Curvature_Higgs_CT_L4[5][5][6][6] = dL2;
    Curvature_Higgs_CT_L4[5][5][7][7] = dL2;
    Curvature_Higgs_CT_L4[5][5][8][8] = dL8 / 0.2e1;
    Curvature_Higgs_CT_L4[5][5][9][9] = dL8 / 0.2e1;
    Curvature_Higgs_CT_L4[5][6][0][2] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[5][6][0][3] = -dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[5][6][1][2] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[5][6][1][3] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[5][6][2][0] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[5][6][2][1] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[5][6][3][0] = -dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[5][6][3][1] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[5][6][5][6] = dL2;
    Curvature_Higgs_CT_L4[5][6][6][5] = dL2;
    Curvature_Higgs_CT_L4[5][7][0][2] = dL4 / 0.2e1 - dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[5][7][0][3] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[5][7][1][2] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[5][7][1][3] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[5][7][2][0] = dL4 / 0.2e1 - dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[5][7][2][1] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[5][7][3][0] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[5][7][3][1] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[5][7][5][7] = dL2;
    Curvature_Higgs_CT_L4[5][7][7][5] = dL2;
    Curvature_Higgs_CT_L4[5][8][5][8] = dL8 / 0.2e1;
    Curvature_Higgs_CT_L4[5][8][8][5] = dL8 / 0.2e1;
    Curvature_Higgs_CT_L4[5][9][5][9] = dL8 / 0.2e1;
    Curvature_Higgs_CT_L4[5][9][9][5] = dL8 / 0.2e1;
    Curvature_Higgs_CT_L4[6][0][0][6] = dL3;
    Curvature_Higgs_CT_L4[6][0][2][4] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[6][0][2][5] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[6][0][3][4] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[6][0][3][5] = -dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[6][0][4][2] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[6][0][4][3] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[6][0][5][2] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[6][0][5][3] = -dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[6][0][6][0] = dL3;
    Curvature_Higgs_CT_L4[6][1][1][6] = dL3;
    Curvature_Higgs_CT_L4[6][1][2][4] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[6][1][2][5] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[6][1][3][4] = dL4 / 0.2e1 - dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[6][1][3][5] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[6][1][4][2] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[6][1][4][3] = dL4 / 0.2e1 - dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[6][1][5][2] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[6][1][5][3] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[6][1][6][1] = dL3;
    Curvature_Higgs_CT_L4[6][2][0][4] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[6][2][0][5] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[6][2][1][4] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[6][2][1][5] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[6][2][2][6] = dL3 + dL4 + dReL5;
    Curvature_Higgs_CT_L4[6][2][2][7] = -dImL5;
    Curvature_Higgs_CT_L4[6][2][3][6] = dImL5;
    Curvature_Higgs_CT_L4[6][2][3][7] = dReL5;
    Curvature_Higgs_CT_L4[6][2][4][0] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[6][2][4][1] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[6][2][5][0] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[6][2][5][1] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[6][2][6][2] = dL3 + dL4 + dReL5;
    Curvature_Higgs_CT_L4[6][2][6][3] = dImL5;
    Curvature_Higgs_CT_L4[6][2][7][2] = -dImL5;
    Curvature_Higgs_CT_L4[6][2][7][3] = dReL5;
    Curvature_Higgs_CT_L4[6][3][0][4] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[6][3][0][5] = -dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[6][3][1][4] = dL4 / 0.2e1 - dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[6][3][1][5] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[6][3][2][6] = dImL5;
    Curvature_Higgs_CT_L4[6][3][2][7] = dReL5;
    Curvature_Higgs_CT_L4[6][3][3][6] = dL3 + dL4 - dReL5;
    Curvature_Higgs_CT_L4[6][3][3][7] = dImL5;
    Curvature_Higgs_CT_L4[6][3][4][0] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[6][3][4][1] = dL4 / 0.2e1 - dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[6][3][5][0] = -dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[6][3][5][1] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[6][3][6][2] = dImL5;
    Curvature_Higgs_CT_L4[6][3][6][3] = dL3 + dL4 - dReL5;
    Curvature_Higgs_CT_L4[6][3][7][2] = dReL5;
    Curvature_Higgs_CT_L4[6][3][7][3] = dImL5;
    Curvature_Higgs_CT_L4[6][4][0][2] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[6][4][0][3] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[6][4][1][2] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[6][4][1][3] = dL4 / 0.2e1 - dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[6][4][2][0] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[6][4][2][1] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[6][4][3][0] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[6][4][3][1] = dL4 / 0.2e1 - dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[6][4][4][6] = dL2;
    Curvature_Higgs_CT_L4[6][4][6][4] = dL2;
    Curvature_Higgs_CT_L4[6][5][0][2] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[6][5][0][3] = -dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[6][5][1][2] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[6][5][1][3] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[6][5][2][0] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[6][5][2][1] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[6][5][3][0] = -dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[6][5][3][1] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[6][5][5][6] = dL2;
    Curvature_Higgs_CT_L4[6][5][6][5] = dL2;
    Curvature_Higgs_CT_L4[6][6][0][0] = dL3;
    Curvature_Higgs_CT_L4[6][6][1][1] = dL3;
    Curvature_Higgs_CT_L4[6][6][2][2] = dL3 + dL4 + dReL5;
    Curvature_Higgs_CT_L4[6][6][2][3] = dImL5;
    Curvature_Higgs_CT_L4[6][6][3][2] = dImL5;
    Curvature_Higgs_CT_L4[6][6][3][3] = dL3 + dL4 - dReL5;
    Curvature_Higgs_CT_L4[6][6][4][4] = dL2;
    Curvature_Higgs_CT_L4[6][6][5][5] = dL2;
    Curvature_Higgs_CT_L4[6][6][6][6] = 3 * dL2;
    Curvature_Higgs_CT_L4[6][6][7][7] = dL2;
    Curvature_Higgs_CT_L4[6][6][8][8] = dL8 / 0.2e1;
    Curvature_Higgs_CT_L4[6][6][9][9] = dL8 / 0.2e1;
    Curvature_Higgs_CT_L4[6][7][2][2] = -dImL5;
    Curvature_Higgs_CT_L4[6][7][2][3] = dReL5;
    Curvature_Higgs_CT_L4[6][7][3][2] = dReL5;
    Curvature_Higgs_CT_L4[6][7][3][3] = dImL5;
    Curvature_Higgs_CT_L4[6][7][6][7] = dL2;
    Curvature_Higgs_CT_L4[6][7][7][6] = dL2;
    Curvature_Higgs_CT_L4[6][8][6][8] = dL8 / 0.2e1;
    Curvature_Higgs_CT_L4[6][8][8][6] = dL8 / 0.2e1;
    Curvature_Higgs_CT_L4[6][9][6][9] = dL8 / 0.2e1;
    Curvature_Higgs_CT_L4[6][9][9][6] = dL8 / 0.2e1;
    Curvature_Higgs_CT_L4[7][0][0][7] = dL3;
    Curvature_Higgs_CT_L4[7][0][2][4] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[7][0][2][5] = dL4 / 0.2e1 - dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[7][0][3][4] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[7][0][3][5] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[7][0][4][2] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[7][0][4][3] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[7][0][5][2] = dL4 / 0.2e1 - dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[7][0][5][3] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[7][0][7][0] = dL3;
    Curvature_Higgs_CT_L4[7][1][1][7] = dL3;
    Curvature_Higgs_CT_L4[7][1][2][4] = -dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[7][1][2][5] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[7][1][3][4] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[7][1][3][5] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[7][1][4][2] = -dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[7][1][4][3] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[7][1][5][2] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[7][1][5][3] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[7][1][7][1] = dL3;
    Curvature_Higgs_CT_L4[7][2][0][4] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[7][2][0][5] = dL4 / 0.2e1 - dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[7][2][1][4] = -dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[7][2][1][5] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[7][2][2][6] = -dImL5;
    Curvature_Higgs_CT_L4[7][2][2][7] = dL3 + dL4 - dReL5;
    Curvature_Higgs_CT_L4[7][2][3][6] = dReL5;
    Curvature_Higgs_CT_L4[7][2][3][7] = -dImL5;
    Curvature_Higgs_CT_L4[7][2][4][0] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[7][2][4][1] = -dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[7][2][5][0] = dL4 / 0.2e1 - dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[7][2][5][1] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[7][2][6][2] = -dImL5;
    Curvature_Higgs_CT_L4[7][2][6][3] = dReL5;
    Curvature_Higgs_CT_L4[7][2][7][2] = dL3 + dL4 - dReL5;
    Curvature_Higgs_CT_L4[7][2][7][3] = -dImL5;
    Curvature_Higgs_CT_L4[7][3][0][4] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[7][3][0][5] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[7][3][1][4] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[7][3][1][5] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[7][3][2][6] = dReL5;
    Curvature_Higgs_CT_L4[7][3][2][7] = -dImL5;
    Curvature_Higgs_CT_L4[7][3][3][6] = dImL5;
    Curvature_Higgs_CT_L4[7][3][3][7] = dL3 + dL4 + dReL5;
    Curvature_Higgs_CT_L4[7][3][4][0] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[7][3][4][1] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[7][3][5][0] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[7][3][5][1] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[7][3][6][2] = dReL5;
    Curvature_Higgs_CT_L4[7][3][6][3] = dImL5;
    Curvature_Higgs_CT_L4[7][3][7][2] = -dImL5;
    Curvature_Higgs_CT_L4[7][3][7][3] = dL3 + dL4 + dReL5;
    Curvature_Higgs_CT_L4[7][4][0][2] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[7][4][0][3] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[7][4][1][2] = -dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[7][4][1][3] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[7][4][2][0] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[7][4][2][1] = -dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[7][4][3][0] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[7][4][3][1] = dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[7][4][4][7] = dL2;
    Curvature_Higgs_CT_L4[7][4][7][4] = dL2;
    Curvature_Higgs_CT_L4[7][5][0][2] = dL4 / 0.2e1 - dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[7][5][0][3] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[7][5][1][2] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[7][5][1][3] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[7][5][2][0] = dL4 / 0.2e1 - dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[7][5][2][1] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[7][5][3][0] = -dImL5 / 0.2e1;
    Curvature_Higgs_CT_L4[7][5][3][1] = dL4 / 0.2e1 + dReL5 / 0.2e1;
    Curvature_Higgs_CT_L4[7][5][5][7] = dL2;
    Curvature_Higgs_CT_L4[7][5][7][5] = dL2;
    Curvature_Higgs_CT_L4[7][6][2][2] = -dImL5;
    Curvature_Higgs_CT_L4[7][6][2][3] = dReL5;
    Curvature_Higgs_CT_L4[7][6][3][2] = dReL5;
    Curvature_Higgs_CT_L4[7][6][3][3] = dImL5;
    Curvature_Higgs_CT_L4[7][6][6][7] = dL2;
    Curvature_Higgs_CT_L4[7][6][7][6] = dL2;
    Curvature_Higgs_CT_L4[7][7][0][0] = dL3;
    Curvature_Higgs_CT_L4[7][7][1][1] = dL3;
    Curvature_Higgs_CT_L4[7][7][2][2] = dL3 + dL4 - dReL5;
    Curvature_Higgs_CT_L4[7][7][2][3] = -dImL5;
    Curvature_Higgs_CT_L4[7][7][3][2] = -dImL5;
    Curvature_Higgs_CT_L4[7][7][3][3] = dL3 + dL4 + dReL5;
    Curvature_Higgs_CT_L4[7][7][4][4] = dL2;
    Curvature_Higgs_CT_L4[7][7][5][5] = dL2;
    Curvature_Higgs_CT_L4[7][7][6][6] = dL2;
    Curvature_Higgs_CT_L4[7][7][7][7] = 3 * dL2;
    Curvature_Higgs_CT_L4[7][7][8][8] = dL8 / 0.2e1;
    Curvature_Higgs_CT_L4[7][7][9][9] = dL8 / 0.2e1;
    Curvature_Higgs_CT_L4[7][8][7][8] = dL8 / 0.2e1;
    Curvature_Higgs_CT_L4[7][8][8][7] = dL8 / 0.2e1;
    Curvature_Higgs_CT_L4[7][9][7][9] = dL8 / 0.2e1;
    Curvature_Higgs_CT_L4[7][9][9][7] = dL8 / 0.2e1;
    Curvature_Higgs_CT_L4[8][0][0][8] = dL7 / 0.2e1;
    Curvature_Higgs_CT_L4[8][0][8][0] = dL7 / 0.2e1;
    Curvature_Higgs_CT_L4[8][1][1][8] = dL7 / 0.2e1;
    Curvature_Higgs_CT_L4[8][1][8][1] = dL7 / 0.2e1;
    Curvature_Higgs_CT_L4[8][2][2][8] = dL7 / 0.2e1;
    Curvature_Higgs_CT_L4[8][2][8][2] = dL7 / 0.2e1;
    Curvature_Higgs_CT_L4[8][3][3][8] = dL7 / 0.2e1;
    Curvature_Higgs_CT_L4[8][3][8][3] = dL7 / 0.2e1;
    Curvature_Higgs_CT_L4[8][4][4][8] = dL8 / 0.2e1;
    Curvature_Higgs_CT_L4[8][4][8][4] = dL8 / 0.2e1;
    Curvature_Higgs_CT_L4[8][5][5][8] = dL8 / 0.2e1;
    Curvature_Higgs_CT_L4[8][5][8][5] = dL8 / 0.2e1;
    Curvature_Higgs_CT_L4[8][6][6][8] = dL8 / 0.2e1;
    Curvature_Higgs_CT_L4[8][6][8][6] = dL8 / 0.2e1;
    Curvature_Higgs_CT_L4[8][7][7][8] = dL8 / 0.2e1;
    Curvature_Higgs_CT_L4[8][7][8][7] = dL8 / 0.2e1;
    Curvature_Higgs_CT_L4[8][8][0][0] = dL7 / 0.2e1;
    Curvature_Higgs_CT_L4[8][8][1][1] = dL7 / 0.2e1;
    Curvature_Higgs_CT_L4[8][8][2][2] = dL7 / 0.2e1;
    Curvature_Higgs_CT_L4[8][8][3][3] = dL7 / 0.2e1;
    Curvature_Higgs_CT_L4[8][8][4][4] = dL8 / 0.2e1;
    Curvature_Higgs_CT_L4[8][8][5][5] = dL8 / 0.2e1;
    Curvature_Higgs_CT_L4[8][8][6][6] = dL8 / 0.2e1;
    Curvature_Higgs_CT_L4[8][8][7][7] = dL8 / 0.2e1;
    Curvature_Higgs_CT_L4[8][8][8][8] = 0.3e1 / 0.4e1 * dL6;
    Curvature_Higgs_CT_L4[8][8][9][9] = dL6 / 0.4e1;
    Curvature_Higgs_CT_L4[8][9][8][9] = dL6 / 0.4e1;
    Curvature_Higgs_CT_L4[8][9][9][8] = dL6 / 0.4e1;
    Curvature_Higgs_CT_L4[9][0][0][9] = dL7 / 0.2e1;
    Curvature_Higgs_CT_L4[9][0][9][0] = dL7 / 0.2e1;
    Curvature_Higgs_CT_L4[9][1][1][9] = dL7 / 0.2e1;
    Curvature_Higgs_CT_L4[9][1][9][1] = dL7 / 0.2e1;
    Curvature_Higgs_CT_L4[9][2][2][9] = dL7 / 0.2e1;
    Curvature_Higgs_CT_L4[9][2][9][2] = dL7 / 0.2e1;
    Curvature_Higgs_CT_L4[9][3][3][9] = dL7 / 0.2e1;
    Curvature_Higgs_CT_L4[9][3][9][3] = dL7 / 0.2e1;
    Curvature_Higgs_CT_L4[9][4][4][9] = dL8 / 0.2e1;
    Curvature_Higgs_CT_L4[9][4][9][4] = dL8 / 0.2e1;
    Curvature_Higgs_CT_L4[9][5][5][9] = dL8 / 0.2e1;
    Curvature_Higgs_CT_L4[9][5][9][5] = dL8 / 0.2e1;
    Curvature_Higgs_CT_L4[9][6][6][9] = dL8 / 0.2e1;
    Curvature_Higgs_CT_L4[9][6][9][6] = dL8 / 0.2e1;
    Curvature_Higgs_CT_L4[9][7][7][9] = dL8 / 0.2e1;
    Curvature_Higgs_CT_L4[9][7][9][7] = dL8 / 0.2e1;
    Curvature_Higgs_CT_L4[9][8][8][9] = dL6 / 0.4e1;
    Curvature_Higgs_CT_L4[9][8][9][8] = dL6 / 0.4e1;
    Curvature_Higgs_CT_L4[9][9][0][0] = dL7 / 0.2e1;
    Curvature_Higgs_CT_L4[9][9][1][1] = dL7 / 0.2e1;
    Curvature_Higgs_CT_L4[9][9][2][2] = dL7 / 0.2e1;
    Curvature_Higgs_CT_L4[9][9][3][3] = dL7 / 0.2e1;
    Curvature_Higgs_CT_L4[9][9][4][4] = dL8 / 0.2e1;
    Curvature_Higgs_CT_L4[9][9][5][5] = dL8 / 0.2e1;
    Curvature_Higgs_CT_L4[9][9][6][6] = dL8 / 0.2e1;
    Curvature_Higgs_CT_L4[9][9][7][7] = dL8 / 0.2e1;
    Curvature_Higgs_CT_L4[9][9][8][8] = dL6 / 0.4e1;
    Curvature_Higgs_CT_L4[9][9][9][9] = 0.3e1 / 0.4e1 * dL6;

    //std::cout << "Curvature tensor 4 which should NOT be 0 = " << Curvature_Higgs_CT_L4[9][9][8][8] << "\n" << std::endl;
    //std::cout << "Curvature tensor 4 which should be 0 = " << Curvature_Higgs_CT_L4[9][9][8][7] << "\n" << std::endl;

    double CurvatureL1HiggsCT;
    double CurvatureL2HiggsCT;
    double CurvatureL4HiggsCT;
    CurvatureL1HiggsCT = 0;
    CurvatureL2HiggsCT = 0;
    CurvatureL4HiggsCT = 0;

    for(int i=0;i<10;i++){
        CurvatureL1HiggsCT += Curvature_Higgs_CT_L1[i];
        for(int j=0;j<10;j++){
            CurvatureL2HiggsCT += Curvature_Higgs_CT_L2[i][j];
            for(int k=0;k<10;k++){
                for(int l=0;l<10;l++){
                    CurvatureL4HiggsCT += Curvature_Higgs_CT_L4[i][j][k][l];
                }
            }
        }
    }

    std::cout << "Curvature tensor Higgs L1 added " << CurvatureL1HiggsCT << "\n" << std::endl;
    std::cout << "Curvature tensor Higgs L2 added " << CurvatureL2HiggsCT << "\n" << std::endl;
    std::cout << "Curvature tensor Higgs L4 added " << CurvatureL4HiggsCT << "\n" << std::endl;



}



/**
 * console output of all Parameters
 */
void Class_Potential_CN2HDM::write() const {

    // Writing out input prameters for user.

	std::cout << "The parameters are : " << std::endl;
	std::cout << "Type = " << Type << std::endl;
	std::cout << "Re(m_{12}^2) = " << Rem122 << " GeV^2 \n"
            << "Im(m_{12}^2) = " << Imm122 << " GeV^2 \n"
            << "m_{11}^2 = " << m112 << " GeV^2 \n"
            << "m_{22}^2 = " << m222 << " GeV^2 \n"
            << "m_s^2 = " << ms2 << " GeV^2 \n"
            << "m_{DM}^2 = " << mds << " GeV^2 \n"
            << "v_s = " << vs << " GeV \n"
            << "v = " << C_vev0 << " GeV \n"
            << "Tan(Beta) = " << TanBeta << "\n"
            << "L1 = " << L1 << "\n"
            << "L2 = " << L2 << "\n"
            << "L3 = " << L3 << "\n"
            << "L4 = " << L4 << "\n"
            << "Re(L5) = " << ReL5 << "\n"
            << "Im(L5) = " << ImL5 << "\n"
            << "L6 = " << L6 << "\n"
            << "L7 = " << L7 << "\n"
            << "L8 = " << L8 << std::endl;

    std::cout << "The CounterTerm parameters are : " << std::endl;
    std::cout << "dRe(m_{12}^2) = " << dRem122 << " GeV^2 \n"
            << "dIm(m_{12}^2) = " << dImm122 << " GeV^2 \n"
            << "dm_{11}^2 = " << dm112 << " GeV^2 \n"
            << "dm_{22}^2 = " << dm222 << " GeV^2 \n"
            << "dm_s^2 = " << dms2 << " GeV^2 \n"
            << "dm_{DM}^2 = " << dmDM2 << " GeV^2 \n"
            << "dL1 = " << dL1 << "\n"
            << "dL2 = " << dL2 << "\n"
            << "dL3 = " << dL3 << "\n"
            << "dL4 = " << dL4 << "\n"
            << "dRe(L5) = " << dReL5 << "\n"
            << "dIm(L5) = " << dImL5 << "\n"
            << "dL6 = " << dL6 << "\n"
            << "dL7 = " << dL7 << "\n"
            << "dL8 = " << dL8 << "\n"
            << "dT1 = " << dT1 << "\n"
            << "dT2 = " << dT2 << "\n"
            << "dT3 = " << dT3 << "\n"
            << "dT4 = " << dT4 << "\n"
            << "dT5 = " << dT5 << "\n"
            << "dT6 = " << dT6 << "\n"
            << "dT7 = " << dT7 << "\n"
            << "dT8 = " << dT8 << "\n"
            << "dT9 = " << dT9 << "\n"
            << "dT10 = " << dT10 << std::endl;

	std::cout << "The scale is given by mu = " << scale << " GeV " << std::endl;

    std::cout << "The masses are given by :" << std::endl;
    std::vector<double> mhiggssquared;

    mhiggssquared=HiggsMassesSquared(vevTree,0,0);

    for(size_t i=0;i<NHiggs;i++){
        double mtmp = mhiggssquared.at(i);
        std::cout << "m_{h_" << i << "}^2 = " << mtmp << "\t m_{h_" << i << "} = ";
        if(mtmp < 0 ) std::cout << " i * ";
        std::cout << std::sqrt(std::abs(mtmp)) << " GeV " << std::endl;
    }



}


/**
 * Calculates the counterterms. Here you need to work out the scheme and implement the formulas.
 */
std::vector<double> Class_Potential_CN2HDM::calc_CT() const{
    using std::pow;
    std::vector<double> parCT;
	bool Debug=false;
	if(Debug) std::cout << "Start " << __func__ << std::endl;

    if(!SetCurvatureDone){
        std::string retmes = __func__;
        retmes += " was called before SetCurvatureArrays()!\n";
        throw std::runtime_error(retmes);
    }
    if(!CalcCouplingsdone){
        std::string retmes = __func__;
        retmes += " was called before CalculatePhysicalCouplings()!\n";
        throw std::runtime_error(retmes);
    }

	if(Debug) {
	std::cout << "Couplings done " << std::endl;
	}
    std::vector<double> WeinbergNabla,WeinbergHesse;
    WeinbergNabla=WeinbergFirstDerivative();
    WeinbergHesse=WeinbergSecondDerivative();

	if(Debug) std::cout << "Finished Nabla/Weinberg Derivatives " << std::endl;
    //std::cout << "1st Nabla Weinberg Derviative " << WeinbergNabla[1] << std::endl;
    //std::cout << "if the same as Phillip's should just give 0 " << std::endl;

	VectorXd NablaWeinberg(NHiggs);
	MatrixXd HesseWeinberg(NHiggs,NHiggs),HiggsRot(NHiggs,NHiggs);
	for(size_t i=0;i<NHiggs;i++)
	{
		NablaWeinberg[i] = WeinbergNabla[i];
        //std::cout << NablaWeinberg[i] << std::endl;
		for(size_t j=0;j<NHiggs;j++) HesseWeinberg(i,j) = WeinbergHesse.at(j*NHiggs+i);
	}

    if(Debug){
		std::cout << "NablaWeinberg = " << NablaWeinberg.transpose() << std::endl;
		std::cout << "Hessian Weinberg : " << std::endl
				<< HesseWeinberg << std::endl;
	}

    std::cout << "Printing NW and HW terms \n" << std::endl;
    std::cout << "HesseWeinberg(0, 0) = " << HesseWeinberg(0, 0) << "\n" << std::endl;
    std::cout << "S HesseWeinberg(3, 3) = P HesseWeinberg(6, 6) " << HesseWeinberg(3, 3) << "\n" << std::endl;
    std::cout << "S HesseWeinberg(7, 7) = P HesseWeinberg(7, 7) " << HesseWeinberg(7, 7) << "\n" << std::endl;
    std::cout << "S HesseWeinberg(2, 6) = P HesseWeinberg(4, 5) " << HesseWeinberg(2, 6) << "\n" << std::endl;
    std::cout << "S HesseWeinberg(3, 7) = P HesseWeinberg(6, 7) " << HesseWeinberg(3, 7) << "\n" << std::endl;
    std::cout << "S HesseWeinberg(6, 6) = P HesseWeinberg(5, 5) " << HesseWeinberg(6, 6) << "\n" << std::endl;
    std::cout << "S HesseWeinberg(6, 8) = P HesseWeinberg(5, 8) " << HesseWeinberg(6, 8) << "\n" << std::endl;

    std::cout << "S HesseWeinberg(1, 1) = P HesseWeinberg(2, 2) " << HesseWeinberg(1, 1) << "\n" << std::endl;
    std::cout << "S HesseWeinberg(2, 2) = P HesseWeinberg(4, 4) " << HesseWeinberg(2, 2) << "\n" << std::endl;
    std::cout << "S HesseWeinberg(3, 3) = P HesseWeinberg(6, 6) " << HesseWeinberg(3, 3) << "\n" << std::endl;
    std::cout << "S HesseWeinberg(4, 4) = P HesseWeinberg(1, 1) " << HesseWeinberg(4, 4) << "\n" << std::endl;
    std::cout << "S HesseWeinberg(5, 5) = P HesseWeinberg(3, 3) " << HesseWeinberg(5, 5) << "\n" << std::endl;
    std::cout << "S HesseWeinberg(6, 6) = P HesseWeinberg(5, 5) " << HesseWeinberg(6, 6) << "\n" << std::endl;
    std::cout << "S HesseWeinberg(8, 8) = P HesseWeinberg(8, 8) " << HesseWeinberg(8, 8) << "\n" << std::endl;
    std::cout << "S HesseWeinberg(9, 9) = P HesseWeinberg(9, 9) " << HesseWeinberg(9, 9) << "\n" << std::endl;

    std::cout << "S NablaWeinberg(0) = P NablaWeinberg(0) " << NablaWeinberg(0) << "\n" << std::endl;
    std::cout << "S NablaWeinberg(1) = P NablaWeinberg(2) " << NablaWeinberg(1) << "\n" << std::endl;
    std::cout << "S NablaWeinberg(2) = P NablaWeinberg(4) " << NablaWeinberg(2) << "\n" << std::endl;
    std::cout << "S NablaWeinberg(3) = P NablaWeinberg(6) " << NablaWeinberg(3) << "\n" << std::endl;
    std::cout << "S NablaWeinberg(4) = P NablaWeinberg(1) " << NablaWeinberg(4) << "\n" << std::endl;
    std::cout << "S NablaWeinberg(5) = P NablaWeinberg(3) " << NablaWeinberg(5) << "\n" << std::endl;
    std::cout << "S NablaWeinberg(6) = P NablaWeinberg(5) " << NablaWeinberg(6) << "\n" << std::endl;
    std::cout << "S NablaWeinberg(7) = P NablaWeinberg(7) " << NablaWeinberg(7) << "\n" << std::endl;
    std::cout << "S NablaWeinberg(8) = P NablaWeinberg(8) " << NablaWeinberg(8) << "\n" << std::endl;
    std::cout << "S NablaWeinberg(9) = P NablaWeinberg(9) " << NablaWeinberg(9) << "\n" << std::endl;


	// Here you have to use your formulas for the counterterm scheme

    //dm112
    parCT.push_back((-4 * v1 * HesseWeinberg(0, 0) + HesseWeinberg(3, 3) * v1 + HesseWeinberg(2, 2) * v1 + vs * HesseWeinberg(2, 8)
                     + HesseWeinberg(2, 6) * v2 - HesseWeinberg(3, 7) * v2) / v1 / 2);
    //dm222
    parCT.push_back((-4 * v1 * v1 * HesseWeinberg(0, 0) + 4 * v1 * v1 * HesseWeinberg(3, 3) - 3 * HesseWeinberg(7, 7) * v2 * v2
                     + HesseWeinberg(2, 6) * v1 * v2 - HesseWeinberg(3, 7) * v1 * v2 + HesseWeinberg(6, 6) * v2 * v2
                     + vs * HesseWeinberg(6, 8) * v2) * pow(v2, -2) / 2);
    //dRem122
    parCT.push_back(-2 * v1 / v2 * (HesseWeinberg(0, 0) - HesseWeinberg(3, 3)) + HesseWeinberg(3, 7));
    //dImm122
    parCT.push_back((-2 * HesseWeinberg(2, 3) * v1 + HesseWeinberg(3, 6) * v2) / v2);
    //dL1
    parCT.push_back((2 * HesseWeinberg(0, 0) - HesseWeinberg(3, 3) - HesseWeinberg(2, 2)) * pow(v1, -2));
    //dL2
    parCT.push_back((2 * v1 * v1 * (HesseWeinberg(0, 0) - HesseWeinberg(3, 3)) * pow(v2, -2) - HesseWeinberg(6, 6) + HesseWeinberg(7, 7)) * pow(v2, -2));
    //dL3
    parCT.push_back((-HesseWeinberg(2, 6) + HesseWeinberg(3, 7)) / v1 / v2);
    //dL4
    parCT.push_back(0);
    //dReL5
    parCT.push_back((-2 * HesseWeinberg(0, 0) + 2 * HesseWeinberg(3, 3)) * pow(v2, -2));
    //dImL5
    parCT.push_back(-2 * HesseWeinberg(2, 3) * pow(v2, -2));
    //dms2
    parCT.push_back((-HesseWeinberg(9, 9) * vs + v1 * HesseWeinberg(2, 8) + v2 * HesseWeinberg(6, 8) + HesseWeinberg(8, 8) * vs - 2 * NablaWeinberg(8)) / vs);
    //dL6
    parCT.push_back((-4 * HesseWeinberg(8, 8) * vs + 4 * NablaWeinberg(8)) * pow(vs, -3));
    //dL7
    parCT.push_back(-2 * HesseWeinberg(2, 8) / v1 / vs);
    //dL8
    parCT.push_back(-2 * HesseWeinberg(6, 8) / v2 / vs);
    //dmDM2
    parCT.push_back((-HesseWeinberg(9, 9) * vs + NablaWeinberg(8)) / vs);
    //dT1
    parCT.push_back(-NablaWeinberg(0));
    //dT2
    parCT.push_back(-NablaWeinberg(1));
    //dT3
    parCT.push_back(HesseWeinberg(3, 3) * v1 + HesseWeinberg(3, 7) * v2 - NablaWeinberg(2));
    //dT4
    parCT.push_back(-HesseWeinberg(2, 3) * v1 + HesseWeinberg(3, 6) * v2 - NablaWeinberg(3));
    //dT5
    parCT.push_back(-NablaWeinberg(4));
    //dT6
    parCT.push_back(-NablaWeinberg(5));
    //dT7
    parCT.push_back(HesseWeinberg(3, 7) * v1 + HesseWeinberg(7, 7) * v2 - NablaWeinberg(6));
    //dT8
    parCT.push_back((HesseWeinberg(2, 3) * v1 * v1 - v2 * (HesseWeinberg(3, 6) * v1 + NablaWeinberg(7))) / v2);
    //dT9
    parCT.push_back(0);
    //dT10
    parCT.push_back(-NablaWeinberg(9));

    return  parCT;
}




void Class_Potential_CN2HDM::TripleHiggsCouplings()
{
    bool Debug=true;

    if(Debug) std::cout << "in triple higgs coupling section " << std::endl;

	if(Debug) std::cout << "Debug turned on in " << __func__ << std::endl;

	if(!SetCurvatureDone)SetCurvatureArrays();
	if(!CalcCouplingsdone)CalculatePhysicalCouplings();

    MatrixXd HiggsRot(NHiggs,NHiggs);
    for(size_t i=0;i<NHiggs;i++)
    {
        for(size_t j=0;j<NHiggs;j++)
        {
            HiggsRot(i,j) = HiggsRotationMatrix[i][j];
        }
    }

    if(Debug) std::cout << "Rotation matrix text 1,1 component " << HiggsRot(1,1) << std::endl;




    int pos_rho1 = 0, pos_eta1=1, pos_xi1 = 2, pos_psi1 = 3;
    int pos_rho2 = 4, pos_eta2=5, pos_xi2 = 6, pos_psi2 =7;
    int pos_res = 8, pos_ims = 9;

    int posG1 =0, posG2 = 0, posG0 = 0, posMHCS1 =0, posMHCS2 = 0, posDM = 0, posh1 = 0, posh2 = 0, posh3 = 0, posh4=0;

    double testsum = 0;


    std::vector<double> posN(4);
	int countposN = 0;

	std::vector<double> mhiggssquared;
    mhiggssquared=HiggsMassesSquared(vevTree,0,0);
	int nMassless = 0;
	for(auto ms: mhiggssquared) nMassless += (ms == 0);

    for(int i=0;i<nMassless;i++){
        testsum = std::abs(HiggsRot(i,pos_rho1)) + std::abs(HiggsRot(i,pos_rho2));
        if(testsum != 0) posG1 = i;
        testsum = std::abs(HiggsRot(i,pos_eta1)) + std::abs(HiggsRot(i,pos_eta2));
        if(testsum != 0) posG2 = i;
        testsum = std::abs(HiggsRot(i,pos_psi1)) + std::abs(HiggsRot(i,pos_psi2));
        if(testsum != 0) posG0 = i;
        testsum = std::abs(HiggsRot(i,pos_ims));
		if(testsum != 0) posDM = i;

    }

    for(size_t i=nMassless;i<NHiggs;i++)
	{
		testsum = std::abs(HiggsRot(i,pos_ims));
		if(testsum != 0) posDM = i;
		testsum = std::abs(HiggsRot(i,pos_rho1)) + std::abs(HiggsRot(i,pos_rho2));
		if(testsum != 0) posMHCS1 = i;
		testsum = std::abs(HiggsRot(i,pos_eta1)) + std::abs(HiggsRot(i,pos_eta2));
		if(testsum != 0) posMHCS2 = i;
		std::vector<int> neutralFields{pos_xi1,pos_xi2,pos_psi1,pos_psi2,pos_res};
		testsum = 0;
		for(auto field: neutralFields){
			testsum += std::abs(HiggsRot(i,field));
		}
		if(testsum != 0)
		{
			posN[countposN] = i;
			countposN++;
		}
	}

    posh1 = posN[0];
	posh2 = posN[1];
	posh3 = posN[2];
	posh4 = posN[3];




	std::vector<double> HiggsOrder(NHiggs);


	// Here you have to set the vector HiggsOrder. By telling e.g. HiggsOrder[0] = 5 you always want your 6th lightest
	// particle to be the first particle in the vector (which has the index 5 because they are sorted by mass)

	// example for keeping the mass order
    // Set explicitly??

    HiggsOrder[0] = posG1;
	HiggsOrder[1] = posG2;
	HiggsOrder[2] = posMHCS1;
	HiggsOrder[3] = posMHCS2;
	HiggsOrder[4] = posG0;
	HiggsOrder[5] = posh1;
	HiggsOrder[6] = posh2;
	HiggsOrder[7] = posh3;
	HiggsOrder[8] = posh4;
    HiggsOrder[9] = posDM;

    // particles[0] = "G^+";
    // particles[1] = "G^-";
    // particles[2] = "H^+";
    // particles[3] = "H^-";
    // particles[4] = "G^0";
    // particles[5] = "h1";
    // particles[6] = "h2";
    // particles[7] = "h3";
    // particles[8] = "h4";
    // particles[9] = "hDM";

    /////

	if(Debug) std::cout << "Calculate Triple Derivative" << std::endl;
	std::vector<double> TripleDeriv;
    TripleDeriv=WeinbergThirdDerivative();
	if(Debug) std::cout << "Finished calculating triple derivatives " << std::endl;
	std::vector<std::vector<std::vector<double>>> GaugeBasis(NHiggs, std::vector<std::vector<double>>(NHiggs,
				std::vector<double>(NHiggs)));
	for(size_t i=0;i<NHiggs;i++)
	  {
		for(size_t j=0;j<NHiggs;j++)
		{
		  for(size_t k=0;k<NHiggs;k++)
			{
			  GaugeBasis[i][j][k] = TripleDeriv.at(i+j*NHiggs+k*NHiggs*NHiggs);
			}
		}
	  }



	MatrixXd HiggsRotSort(NHiggs,NHiggs);

	for(size_t i=0;i<NHiggs;i++)
	{
		HiggsRotSort.row(i) = HiggsRot.row(HiggsOrder[i]);
	}

	TripleHiggsCorrectionsCWPhysical.resize(NHiggs);
	TripleHiggsCorrectionsTreePhysical.resize(NHiggs);
	TripleHiggsCorrectionsCTPhysical.resize(NHiggs);
	for(size_t i=0;i<NHiggs;i++) {
		TripleHiggsCorrectionsTreePhysical[i].resize(NHiggs);
		TripleHiggsCorrectionsCWPhysical[i].resize(NHiggs);
		TripleHiggsCorrectionsCTPhysical[i].resize(NHiggs);
		for(size_t j=0;j<NHiggs;j++) {
			TripleHiggsCorrectionsCWPhysical[i][j].resize(NHiggs);
			TripleHiggsCorrectionsTreePhysical[i][j].resize(NHiggs);
			TripleHiggsCorrectionsCTPhysical[i][j].resize(NHiggs);
		}
	}

    if(Debug) std::cout << "Setup done " << std::endl;


	for(size_t i=0;i<NHiggs;i++)
	  {
		for(size_t j=0;j<NHiggs;j++)
		{
			for(size_t k=0;k<NHiggs;k++)
			{
			  TripleHiggsCorrectionsCWPhysical[i][j][k] = 0;
			  TripleHiggsCorrectionsTreePhysical[i][j][k] = 0;
			  TripleHiggsCorrectionsCTPhysical[i][j][k] = 0;
			  for(size_t l=0;l<NHiggs;l++)
			  {
				  for(size_t m=0;m<NHiggs;m++)
				  {
					  for(size_t n=0;n<NHiggs;n++)
					  {
						  double RotFac = HiggsRotSort(i,l)*HiggsRotSort(j,m)*HiggsRotSort(k,n);
						  TripleHiggsCorrectionsCWPhysical[i][j][k] += RotFac*GaugeBasis[l][m][n];
						  TripleHiggsCorrectionsTreePhysical[i][j][k] += RotFac*LambdaHiggs_3[l][m][n];
						  TripleHiggsCorrectionsCTPhysical[i][j][k] += RotFac*LambdaHiggs_3_CT[l][m][n];

					  }
				  }
			  }
			}
		}
	  }



}

void Class_Potential_CN2HDM::SetCurvatureArrays(){
  /*
   *  Here you have to set the vectors
   *  Curvature_Higgs_L1,Curvature_Higgs_L2,Curvature_Higgs_L3,Curvature_Higgs_L4
   "
   *  Curvature_Gauge_G2H2
   *  Curvature_Quark_F2H1, Curvature_Lepton_F2H1
   *  as described in the potential in the paper.
   */

	initVectors();
	SetCurvatureDone=true;
	for(size_t i=0;i<NHiggs;i++) HiggsVev[i] = vevTree[i];

    //Curvature_Higgs_L1 - NONE
    //Curvature_Higgs_L2

    Curvature_Higgs_L2[0][0] = m112;
    Curvature_Higgs_L2[0][4] = -Rem122;
    Curvature_Higgs_L2[0][5] = Imm122;
    Curvature_Higgs_L2[1][1] = m112;
    Curvature_Higgs_L2[1][4] = -Imm122;
    Curvature_Higgs_L2[1][5] = -Rem122;
    Curvature_Higgs_L2[2][2] = m112;
    Curvature_Higgs_L2[2][6] = -Rem122;
    Curvature_Higgs_L2[2][7] = Imm122;
    Curvature_Higgs_L2[3][3] = m112;
    Curvature_Higgs_L2[3][6] = -Imm122;
    Curvature_Higgs_L2[3][7] = -Rem122;
    Curvature_Higgs_L2[4][0] = -Rem122;
    Curvature_Higgs_L2[4][1] = -Imm122;
    Curvature_Higgs_L2[4][4] = m222;
    Curvature_Higgs_L2[5][0] = Imm122;
    Curvature_Higgs_L2[5][1] = -Rem122;
    Curvature_Higgs_L2[5][5] = m222;
    Curvature_Higgs_L2[6][2] = -Rem122;
    Curvature_Higgs_L2[6][3] = -Imm122;
    Curvature_Higgs_L2[6][6] = m222;
    Curvature_Higgs_L2[7][2] = Imm122;
    Curvature_Higgs_L2[7][3] = -Rem122;
    Curvature_Higgs_L2[7][7] = m222;
    Curvature_Higgs_L2[8][8] = ms2 / 0.2e1 - mDM2 / 0.2e1;
    Curvature_Higgs_L2[9][9] = ms2 / 0.2e1 + mDM2 / 0.2e1;

    //Curvature_Higgs_L3 - None
    //Curvature_Higgs_L4

    Curvature_Higgs_L4[0][0][0][0] = 3 * L1;
    Curvature_Higgs_L4[0][0][1][1] = L1;
    Curvature_Higgs_L4[0][0][2][2] = L1;
    Curvature_Higgs_L4[0][0][3][3] = L1;
    Curvature_Higgs_L4[0][0][4][4] = L3 + L4 + ReL5;
    Curvature_Higgs_L4[0][0][4][5] = -ImL5;
    Curvature_Higgs_L4[0][0][5][4] = -ImL5;
    Curvature_Higgs_L4[0][0][5][5] = L3 + L4 - ReL5;
    Curvature_Higgs_L4[0][0][6][6] = L3;
    Curvature_Higgs_L4[0][0][7][7] = L3;
    Curvature_Higgs_L4[0][0][8][8] = L7 / 0.2e1;
    Curvature_Higgs_L4[0][0][9][9] = L7 / 0.2e1;
    Curvature_Higgs_L4[0][1][0][1] = L1;
    Curvature_Higgs_L4[0][1][1][0] = L1;
    Curvature_Higgs_L4[0][1][4][4] = ImL5;
    Curvature_Higgs_L4[0][1][4][5] = ReL5;
    Curvature_Higgs_L4[0][1][5][4] = ReL5;
    Curvature_Higgs_L4[0][1][5][5] = -ImL5;
    Curvature_Higgs_L4[0][2][0][2] = L1;
    Curvature_Higgs_L4[0][2][2][0] = L1;
    Curvature_Higgs_L4[0][2][4][6] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[0][2][4][7] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[0][2][5][6] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[0][2][5][7] = L4 / 0.2e1 - ReL5 / 0.2e1;
    Curvature_Higgs_L4[0][2][6][4] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[0][2][6][5] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[0][2][7][4] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[0][2][7][5] = L4 / 0.2e1 - ReL5 / 0.2e1;
    Curvature_Higgs_L4[0][3][0][3] = L1;
    Curvature_Higgs_L4[0][3][3][0] = L1;
    Curvature_Higgs_L4[0][3][4][6] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[0][3][4][7] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[0][3][5][6] = -L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[0][3][5][7] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[0][3][6][4] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[0][3][6][5] = -L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[0][3][7][4] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[0][3][7][5] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[0][4][0][4] = L3 + L4 + ReL5;
    Curvature_Higgs_L4[0][4][0][5] = -ImL5;
    Curvature_Higgs_L4[0][4][1][4] = ImL5;
    Curvature_Higgs_L4[0][4][1][5] = ReL5;
    Curvature_Higgs_L4[0][4][2][6] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[0][4][2][7] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[0][4][3][6] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[0][4][3][7] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[0][4][4][0] = L3 + L4 + ReL5;
    Curvature_Higgs_L4[0][4][4][1] = ImL5;
    Curvature_Higgs_L4[0][4][5][0] = -ImL5;
    Curvature_Higgs_L4[0][4][5][1] = ReL5;
    Curvature_Higgs_L4[0][4][6][2] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[0][4][6][3] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[0][4][7][2] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[0][4][7][3] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[0][5][0][4] = -ImL5;
    Curvature_Higgs_L4[0][5][0][5] = L3 + L4 - ReL5;
    Curvature_Higgs_L4[0][5][1][4] = ReL5;
    Curvature_Higgs_L4[0][5][1][5] = -ImL5;
    Curvature_Higgs_L4[0][5][2][6] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[0][5][2][7] = L4 / 0.2e1 - ReL5 / 0.2e1;
    Curvature_Higgs_L4[0][5][3][6] = -L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[0][5][3][7] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[0][5][4][0] = -ImL5;
    Curvature_Higgs_L4[0][5][4][1] = ReL5;
    Curvature_Higgs_L4[0][5][5][0] = L3 + L4 - ReL5;
    Curvature_Higgs_L4[0][5][5][1] = -ImL5;
    Curvature_Higgs_L4[0][5][6][2] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[0][5][6][3] = -L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[0][5][7][2] = L4 / 0.2e1 - ReL5 / 0.2e1;
    Curvature_Higgs_L4[0][5][7][3] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[0][6][0][6] = L3;
    Curvature_Higgs_L4[0][6][2][4] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[0][6][2][5] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[0][6][3][4] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[0][6][3][5] = -L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[0][6][4][2] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[0][6][4][3] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[0][6][5][2] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[0][6][5][3] = -L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[0][6][6][0] = L3;
    Curvature_Higgs_L4[0][7][0][7] = L3;
    Curvature_Higgs_L4[0][7][2][4] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[0][7][2][5] = L4 / 0.2e1 - ReL5 / 0.2e1;
    Curvature_Higgs_L4[0][7][3][4] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[0][7][3][5] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[0][7][4][2] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[0][7][4][3] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[0][7][5][2] = L4 / 0.2e1 - ReL5 / 0.2e1;
    Curvature_Higgs_L4[0][7][5][3] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[0][7][7][0] = L3;
    Curvature_Higgs_L4[0][8][0][8] = L7 / 0.2e1;
    Curvature_Higgs_L4[0][8][8][0] = L7 / 0.2e1;
    Curvature_Higgs_L4[0][9][0][9] = L7 / 0.2e1;
    Curvature_Higgs_L4[0][9][9][0] = L7 / 0.2e1;
    Curvature_Higgs_L4[1][0][0][1] = L1;
    Curvature_Higgs_L4[1][0][1][0] = L1;
    Curvature_Higgs_L4[1][0][4][4] = ImL5;
    Curvature_Higgs_L4[1][0][4][5] = ReL5;
    Curvature_Higgs_L4[1][0][5][4] = ReL5;
    Curvature_Higgs_L4[1][0][5][5] = -ImL5;
    Curvature_Higgs_L4[1][1][0][0] = L1;
    Curvature_Higgs_L4[1][1][1][1] = 3 * L1;
    Curvature_Higgs_L4[1][1][2][2] = L1;
    Curvature_Higgs_L4[1][1][3][3] = L1;
    Curvature_Higgs_L4[1][1][4][4] = L3 + L4 - ReL5;
    Curvature_Higgs_L4[1][1][4][5] = ImL5;
    Curvature_Higgs_L4[1][1][5][4] = ImL5;
    Curvature_Higgs_L4[1][1][5][5] = L3 + L4 + ReL5;
    Curvature_Higgs_L4[1][1][6][6] = L3;
    Curvature_Higgs_L4[1][1][7][7] = L3;
    Curvature_Higgs_L4[1][1][8][8] = L7 / 0.2e1;
    Curvature_Higgs_L4[1][1][9][9] = L7 / 0.2e1;
    Curvature_Higgs_L4[1][2][1][2] = L1;
    Curvature_Higgs_L4[1][2][2][1] = L1;
    Curvature_Higgs_L4[1][2][4][6] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[1][2][4][7] = -L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[1][2][5][6] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[1][2][5][7] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[1][2][6][4] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[1][2][6][5] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[1][2][7][4] = -L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[1][2][7][5] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[1][3][1][3] = L1;
    Curvature_Higgs_L4[1][3][3][1] = L1;
    Curvature_Higgs_L4[1][3][4][6] = L4 / 0.2e1 - ReL5 / 0.2e1;
    Curvature_Higgs_L4[1][3][4][7] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[1][3][5][6] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[1][3][5][7] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[1][3][6][4] = L4 / 0.2e1 - ReL5 / 0.2e1;
    Curvature_Higgs_L4[1][3][6][5] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[1][3][7][4] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[1][3][7][5] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[1][4][0][4] = ImL5;
    Curvature_Higgs_L4[1][4][0][5] = ReL5;
    Curvature_Higgs_L4[1][4][1][4] = L3 + L4 - ReL5;
    Curvature_Higgs_L4[1][4][1][5] = ImL5;
    Curvature_Higgs_L4[1][4][2][6] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[1][4][2][7] = -L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[1][4][3][6] = L4 / 0.2e1 - ReL5 / 0.2e1;
    Curvature_Higgs_L4[1][4][3][7] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[1][4][4][0] = ImL5;
    Curvature_Higgs_L4[1][4][4][1] = L3 + L4 - ReL5;
    Curvature_Higgs_L4[1][4][5][0] = ReL5;
    Curvature_Higgs_L4[1][4][5][1] = ImL5;
    Curvature_Higgs_L4[1][4][6][2] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[1][4][6][3] = L4 / 0.2e1 - ReL5 / 0.2e1;
    Curvature_Higgs_L4[1][4][7][2] = -L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[1][4][7][3] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[1][5][0][4] = ReL5;
    Curvature_Higgs_L4[1][5][0][5] = -ImL5;
    Curvature_Higgs_L4[1][5][1][4] = ImL5;
    Curvature_Higgs_L4[1][5][1][5] = L3 + L4 + ReL5;
    Curvature_Higgs_L4[1][5][2][6] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[1][5][2][7] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[1][5][3][6] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[1][5][3][7] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[1][5][4][0] = ReL5;
    Curvature_Higgs_L4[1][5][4][1] = ImL5;
    Curvature_Higgs_L4[1][5][5][0] = -ImL5;
    Curvature_Higgs_L4[1][5][5][1] = L3 + L4 + ReL5;
    Curvature_Higgs_L4[1][5][6][2] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[1][5][6][3] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[1][5][7][2] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[1][5][7][3] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[1][6][1][6] = L3;
    Curvature_Higgs_L4[1][6][2][4] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[1][6][2][5] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[1][6][3][4] = L4 / 0.2e1 - ReL5 / 0.2e1;
    Curvature_Higgs_L4[1][6][3][5] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[1][6][4][2] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[1][6][4][3] = L4 / 0.2e1 - ReL5 / 0.2e1;
    Curvature_Higgs_L4[1][6][5][2] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[1][6][5][3] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[1][6][6][1] = L3;
    Curvature_Higgs_L4[1][7][1][7] = L3;
    Curvature_Higgs_L4[1][7][2][4] = -L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[1][7][2][5] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[1][7][3][4] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[1][7][3][5] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[1][7][4][2] = -L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[1][7][4][3] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[1][7][5][2] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[1][7][5][3] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[1][7][7][1] = L3;
    Curvature_Higgs_L4[1][8][1][8] = L7 / 0.2e1;
    Curvature_Higgs_L4[1][8][8][1] = L7 / 0.2e1;
    Curvature_Higgs_L4[1][9][1][9] = L7 / 0.2e1;
    Curvature_Higgs_L4[1][9][9][1] = L7 / 0.2e1;
    Curvature_Higgs_L4[2][0][0][2] = L1;
    Curvature_Higgs_L4[2][0][2][0] = L1;
    Curvature_Higgs_L4[2][0][4][6] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[2][0][4][7] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[2][0][5][6] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[2][0][5][7] = L4 / 0.2e1 - ReL5 / 0.2e1;
    Curvature_Higgs_L4[2][0][6][4] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[2][0][6][5] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[2][0][7][4] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[2][0][7][5] = L4 / 0.2e1 - ReL5 / 0.2e1;
    Curvature_Higgs_L4[2][1][1][2] = L1;
    Curvature_Higgs_L4[2][1][2][1] = L1;
    Curvature_Higgs_L4[2][1][4][6] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[2][1][4][7] = -L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[2][1][5][6] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[2][1][5][7] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[2][1][6][4] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[2][1][6][5] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[2][1][7][4] = -L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[2][1][7][5] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[2][2][0][0] = L1;
    Curvature_Higgs_L4[2][2][1][1] = L1;
    Curvature_Higgs_L4[2][2][2][2] = 3 * L1;
    Curvature_Higgs_L4[2][2][3][3] = L1;
    Curvature_Higgs_L4[2][2][4][4] = L3;
    Curvature_Higgs_L4[2][2][5][5] = L3;
    Curvature_Higgs_L4[2][2][6][6] = L3 + L4 + ReL5;
    Curvature_Higgs_L4[2][2][6][7] = -ImL5;
    Curvature_Higgs_L4[2][2][7][6] = -ImL5;
    Curvature_Higgs_L4[2][2][7][7] = L3 + L4 - ReL5;
    Curvature_Higgs_L4[2][2][8][8] = L7 / 0.2e1;
    Curvature_Higgs_L4[2][2][9][9] = L7 / 0.2e1;
    Curvature_Higgs_L4[2][3][2][3] = L1;
    Curvature_Higgs_L4[2][3][3][2] = L1;
    Curvature_Higgs_L4[2][3][6][6] = ImL5;
    Curvature_Higgs_L4[2][3][6][7] = ReL5;
    Curvature_Higgs_L4[2][3][7][6] = ReL5;
    Curvature_Higgs_L4[2][3][7][7] = -ImL5;
    Curvature_Higgs_L4[2][4][0][6] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[2][4][0][7] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[2][4][1][6] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[2][4][1][7] = -L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[2][4][2][4] = L3;
    Curvature_Higgs_L4[2][4][4][2] = L3;
    Curvature_Higgs_L4[2][4][6][0] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[2][4][6][1] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[2][4][7][0] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[2][4][7][1] = -L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[2][5][0][6] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[2][5][0][7] = L4 / 0.2e1 - ReL5 / 0.2e1;
    Curvature_Higgs_L4[2][5][1][6] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[2][5][1][7] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[2][5][2][5] = L3;
    Curvature_Higgs_L4[2][5][5][2] = L3;
    Curvature_Higgs_L4[2][5][6][0] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[2][5][6][1] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[2][5][7][0] = L4 / 0.2e1 - ReL5 / 0.2e1;
    Curvature_Higgs_L4[2][5][7][1] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[2][6][0][4] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[2][6][0][5] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[2][6][1][4] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[2][6][1][5] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[2][6][2][6] = L3 + L4 + ReL5;
    Curvature_Higgs_L4[2][6][2][7] = -ImL5;
    Curvature_Higgs_L4[2][6][3][6] = ImL5;
    Curvature_Higgs_L4[2][6][3][7] = ReL5;
    Curvature_Higgs_L4[2][6][4][0] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[2][6][4][1] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[2][6][5][0] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[2][6][5][1] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[2][6][6][2] = L3 + L4 + ReL5;
    Curvature_Higgs_L4[2][6][6][3] = ImL5;
    Curvature_Higgs_L4[2][6][7][2] = -ImL5;
    Curvature_Higgs_L4[2][6][7][3] = ReL5;
    Curvature_Higgs_L4[2][7][0][4] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[2][7][0][5] = L4 / 0.2e1 - ReL5 / 0.2e1;
    Curvature_Higgs_L4[2][7][1][4] = -L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[2][7][1][5] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[2][7][2][6] = -ImL5;
    Curvature_Higgs_L4[2][7][2][7] = L3 + L4 - ReL5;
    Curvature_Higgs_L4[2][7][3][6] = ReL5;
    Curvature_Higgs_L4[2][7][3][7] = -ImL5;
    Curvature_Higgs_L4[2][7][4][0] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[2][7][4][1] = -L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[2][7][5][0] = L4 / 0.2e1 - ReL5 / 0.2e1;
    Curvature_Higgs_L4[2][7][5][1] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[2][7][6][2] = -ImL5;
    Curvature_Higgs_L4[2][7][6][3] = ReL5;
    Curvature_Higgs_L4[2][7][7][2] = L3 + L4 - ReL5;
    Curvature_Higgs_L4[2][7][7][3] = -ImL5;
    Curvature_Higgs_L4[2][8][2][8] = L7 / 0.2e1;
    Curvature_Higgs_L4[2][8][8][2] = L7 / 0.2e1;
    Curvature_Higgs_L4[2][9][2][9] = L7 / 0.2e1;
    Curvature_Higgs_L4[2][9][9][2] = L7 / 0.2e1;
    Curvature_Higgs_L4[3][0][0][3] = L1;
    Curvature_Higgs_L4[3][0][3][0] = L1;
    Curvature_Higgs_L4[3][0][4][6] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[3][0][4][7] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[3][0][5][6] = -L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[3][0][5][7] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[3][0][6][4] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[3][0][6][5] = -L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[3][0][7][4] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[3][0][7][5] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[3][1][1][3] = L1;
    Curvature_Higgs_L4[3][1][3][1] = L1;
    Curvature_Higgs_L4[3][1][4][6] = L4 / 0.2e1 - ReL5 / 0.2e1;
    Curvature_Higgs_L4[3][1][4][7] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[3][1][5][6] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[3][1][5][7] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[3][1][6][4] = L4 / 0.2e1 - ReL5 / 0.2e1;
    Curvature_Higgs_L4[3][1][6][5] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[3][1][7][4] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[3][1][7][5] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[3][2][2][3] = L1;
    Curvature_Higgs_L4[3][2][3][2] = L1;
    Curvature_Higgs_L4[3][2][6][6] = ImL5;
    Curvature_Higgs_L4[3][2][6][7] = ReL5;
    Curvature_Higgs_L4[3][2][7][6] = ReL5;
    Curvature_Higgs_L4[3][2][7][7] = -ImL5;
    Curvature_Higgs_L4[3][3][0][0] = L1;
    Curvature_Higgs_L4[3][3][1][1] = L1;
    Curvature_Higgs_L4[3][3][2][2] = L1;
    Curvature_Higgs_L4[3][3][3][3] = 3 * L1;
    Curvature_Higgs_L4[3][3][4][4] = L3;
    Curvature_Higgs_L4[3][3][5][5] = L3;
    Curvature_Higgs_L4[3][3][6][6] = L3 + L4 - ReL5;
    Curvature_Higgs_L4[3][3][6][7] = ImL5;
    Curvature_Higgs_L4[3][3][7][6] = ImL5;
    Curvature_Higgs_L4[3][3][7][7] = L3 + L4 + ReL5;
    Curvature_Higgs_L4[3][3][8][8] = L7 / 0.2e1;
    Curvature_Higgs_L4[3][3][9][9] = L7 / 0.2e1;
    Curvature_Higgs_L4[3][4][0][6] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[3][4][0][7] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[3][4][1][6] = L4 / 0.2e1 - ReL5 / 0.2e1;
    Curvature_Higgs_L4[3][4][1][7] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[3][4][3][4] = L3;
    Curvature_Higgs_L4[3][4][4][3] = L3;
    Curvature_Higgs_L4[3][4][6][0] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[3][4][6][1] = L4 / 0.2e1 - ReL5 / 0.2e1;
    Curvature_Higgs_L4[3][4][7][0] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[3][4][7][1] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[3][5][0][6] = -L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[3][5][0][7] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[3][5][1][6] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[3][5][1][7] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[3][5][3][5] = L3;
    Curvature_Higgs_L4[3][5][5][3] = L3;
    Curvature_Higgs_L4[3][5][6][0] = -L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[3][5][6][1] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[3][5][7][0] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[3][5][7][1] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[3][6][0][4] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[3][6][0][5] = -L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[3][6][1][4] = L4 / 0.2e1 - ReL5 / 0.2e1;
    Curvature_Higgs_L4[3][6][1][5] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[3][6][2][6] = ImL5;
    Curvature_Higgs_L4[3][6][2][7] = ReL5;
    Curvature_Higgs_L4[3][6][3][6] = L3 + L4 - ReL5;
    Curvature_Higgs_L4[3][6][3][7] = ImL5;
    Curvature_Higgs_L4[3][6][4][0] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[3][6][4][1] = L4 / 0.2e1 - ReL5 / 0.2e1;
    Curvature_Higgs_L4[3][6][5][0] = -L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[3][6][5][1] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[3][6][6][2] = ImL5;
    Curvature_Higgs_L4[3][6][6][3] = L3 + L4 - ReL5;
    Curvature_Higgs_L4[3][6][7][2] = ReL5;
    Curvature_Higgs_L4[3][6][7][3] = ImL5;
    Curvature_Higgs_L4[3][7][0][4] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[3][7][0][5] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[3][7][1][4] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[3][7][1][5] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[3][7][2][6] = ReL5;
    Curvature_Higgs_L4[3][7][2][7] = -ImL5;
    Curvature_Higgs_L4[3][7][3][6] = ImL5;
    Curvature_Higgs_L4[3][7][3][7] = L3 + L4 + ReL5;
    Curvature_Higgs_L4[3][7][4][0] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[3][7][4][1] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[3][7][5][0] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[3][7][5][1] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[3][7][6][2] = ReL5;
    Curvature_Higgs_L4[3][7][6][3] = ImL5;
    Curvature_Higgs_L4[3][7][7][2] = -ImL5;
    Curvature_Higgs_L4[3][7][7][3] = L3 + L4 + ReL5;
    Curvature_Higgs_L4[3][8][3][8] = L7 / 0.2e1;
    Curvature_Higgs_L4[3][8][8][3] = L7 / 0.2e1;
    Curvature_Higgs_L4[3][9][3][9] = L7 / 0.2e1;
    Curvature_Higgs_L4[3][9][9][3] = L7 / 0.2e1;
    Curvature_Higgs_L4[4][0][0][4] = L3 + L4 + ReL5;
    Curvature_Higgs_L4[4][0][0][5] = -ImL5;
    Curvature_Higgs_L4[4][0][1][4] = ImL5;
    Curvature_Higgs_L4[4][0][1][5] = ReL5;
    Curvature_Higgs_L4[4][0][2][6] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[4][0][2][7] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[4][0][3][6] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[4][0][3][7] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[4][0][4][0] = L3 + L4 + ReL5;
    Curvature_Higgs_L4[4][0][4][1] = ImL5;
    Curvature_Higgs_L4[4][0][5][0] = -ImL5;
    Curvature_Higgs_L4[4][0][5][1] = ReL5;
    Curvature_Higgs_L4[4][0][6][2] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[4][0][6][3] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[4][0][7][2] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[4][0][7][3] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[4][1][0][4] = ImL5;
    Curvature_Higgs_L4[4][1][0][5] = ReL5;
    Curvature_Higgs_L4[4][1][1][4] = L3 + L4 - ReL5;
    Curvature_Higgs_L4[4][1][1][5] = ImL5;
    Curvature_Higgs_L4[4][1][2][6] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[4][1][2][7] = -L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[4][1][3][6] = L4 / 0.2e1 - ReL5 / 0.2e1;
    Curvature_Higgs_L4[4][1][3][7] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[4][1][4][0] = ImL5;
    Curvature_Higgs_L4[4][1][4][1] = L3 + L4 - ReL5;
    Curvature_Higgs_L4[4][1][5][0] = ReL5;
    Curvature_Higgs_L4[4][1][5][1] = ImL5;
    Curvature_Higgs_L4[4][1][6][2] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[4][1][6][3] = L4 / 0.2e1 - ReL5 / 0.2e1;
    Curvature_Higgs_L4[4][1][7][2] = -L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[4][1][7][3] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[4][2][0][6] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[4][2][0][7] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[4][2][1][6] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[4][2][1][7] = -L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[4][2][2][4] = L3;
    Curvature_Higgs_L4[4][2][4][2] = L3;
    Curvature_Higgs_L4[4][2][6][0] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[4][2][6][1] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[4][2][7][0] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[4][2][7][1] = -L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[4][3][0][6] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[4][3][0][7] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[4][3][1][6] = L4 / 0.2e1 - ReL5 / 0.2e1;
    Curvature_Higgs_L4[4][3][1][7] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[4][3][3][4] = L3;
    Curvature_Higgs_L4[4][3][4][3] = L3;
    Curvature_Higgs_L4[4][3][6][0] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[4][3][6][1] = L4 / 0.2e1 - ReL5 / 0.2e1;
    Curvature_Higgs_L4[4][3][7][0] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[4][3][7][1] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[4][4][0][0] = L3 + L4 + ReL5;
    Curvature_Higgs_L4[4][4][0][1] = ImL5;
    Curvature_Higgs_L4[4][4][1][0] = ImL5;
    Curvature_Higgs_L4[4][4][1][1] = L3 + L4 - ReL5;
    Curvature_Higgs_L4[4][4][2][2] = L3;
    Curvature_Higgs_L4[4][4][3][3] = L3;
    Curvature_Higgs_L4[4][4][4][4] = 3 * L2;
    Curvature_Higgs_L4[4][4][5][5] = L2;
    Curvature_Higgs_L4[4][4][6][6] = L2;
    Curvature_Higgs_L4[4][4][7][7] = L2;
    Curvature_Higgs_L4[4][4][8][8] = L8 / 0.2e1;
    Curvature_Higgs_L4[4][4][9][9] = L8 / 0.2e1;
    Curvature_Higgs_L4[4][5][0][0] = -ImL5;
    Curvature_Higgs_L4[4][5][0][1] = ReL5;
    Curvature_Higgs_L4[4][5][1][0] = ReL5;
    Curvature_Higgs_L4[4][5][1][1] = ImL5;
    Curvature_Higgs_L4[4][5][4][5] = L2;
    Curvature_Higgs_L4[4][5][5][4] = L2;
    Curvature_Higgs_L4[4][6][0][2] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[4][6][0][3] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[4][6][1][2] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[4][6][1][3] = L4 / 0.2e1 - ReL5 / 0.2e1;
    Curvature_Higgs_L4[4][6][2][0] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[4][6][2][1] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[4][6][3][0] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[4][6][3][1] = L4 / 0.2e1 - ReL5 / 0.2e1;
    Curvature_Higgs_L4[4][6][4][6] = L2;
    Curvature_Higgs_L4[4][6][6][4] = L2;
    Curvature_Higgs_L4[4][7][0][2] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[4][7][0][3] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[4][7][1][2] = -L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[4][7][1][3] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[4][7][2][0] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[4][7][2][1] = -L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[4][7][3][0] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[4][7][3][1] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[4][7][4][7] = L2;
    Curvature_Higgs_L4[4][7][7][4] = L2;
    Curvature_Higgs_L4[4][8][4][8] = L8 / 0.2e1;
    Curvature_Higgs_L4[4][8][8][4] = L8 / 0.2e1;
    Curvature_Higgs_L4[4][9][4][9] = L8 / 0.2e1;
    Curvature_Higgs_L4[4][9][9][4] = L8 / 0.2e1;
    Curvature_Higgs_L4[5][0][0][4] = -ImL5;
    Curvature_Higgs_L4[5][0][0][5] = L3 + L4 - ReL5;
    Curvature_Higgs_L4[5][0][1][4] = ReL5;
    Curvature_Higgs_L4[5][0][1][5] = -ImL5;
    Curvature_Higgs_L4[5][0][2][6] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[5][0][2][7] = L4 / 0.2e1 - ReL5 / 0.2e1;
    Curvature_Higgs_L4[5][0][3][6] = -L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[5][0][3][7] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[5][0][4][0] = -ImL5;
    Curvature_Higgs_L4[5][0][4][1] = ReL5;
    Curvature_Higgs_L4[5][0][5][0] = L3 + L4 - ReL5;
    Curvature_Higgs_L4[5][0][5][1] = -ImL5;
    Curvature_Higgs_L4[5][0][6][2] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[5][0][6][3] = -L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[5][0][7][2] = L4 / 0.2e1 - ReL5 / 0.2e1;
    Curvature_Higgs_L4[5][0][7][3] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[5][1][0][4] = ReL5;
    Curvature_Higgs_L4[5][1][0][5] = -ImL5;
    Curvature_Higgs_L4[5][1][1][4] = ImL5;
    Curvature_Higgs_L4[5][1][1][5] = L3 + L4 + ReL5;
    Curvature_Higgs_L4[5][1][2][6] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[5][1][2][7] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[5][1][3][6] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[5][1][3][7] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[5][1][4][0] = ReL5;
    Curvature_Higgs_L4[5][1][4][1] = ImL5;
    Curvature_Higgs_L4[5][1][5][0] = -ImL5;
    Curvature_Higgs_L4[5][1][5][1] = L3 + L4 + ReL5;
    Curvature_Higgs_L4[5][1][6][2] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[5][1][6][3] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[5][1][7][2] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[5][1][7][3] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[5][2][0][6] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[5][2][0][7] = L4 / 0.2e1 - ReL5 / 0.2e1;
    Curvature_Higgs_L4[5][2][1][6] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[5][2][1][7] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[5][2][2][5] = L3;
    Curvature_Higgs_L4[5][2][5][2] = L3;
    Curvature_Higgs_L4[5][2][6][0] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[5][2][6][1] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[5][2][7][0] = L4 / 0.2e1 - ReL5 / 0.2e1;
    Curvature_Higgs_L4[5][2][7][1] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[5][3][0][6] = -L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[5][3][0][7] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[5][3][1][6] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[5][3][1][7] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[5][3][3][5] = L3;
    Curvature_Higgs_L4[5][3][5][3] = L3;
    Curvature_Higgs_L4[5][3][6][0] = -L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[5][3][6][1] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[5][3][7][0] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[5][3][7][1] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[5][4][0][0] = -ImL5;
    Curvature_Higgs_L4[5][4][0][1] = ReL5;
    Curvature_Higgs_L4[5][4][1][0] = ReL5;
    Curvature_Higgs_L4[5][4][1][1] = ImL5;
    Curvature_Higgs_L4[5][4][4][5] = L2;
    Curvature_Higgs_L4[5][4][5][4] = L2;
    Curvature_Higgs_L4[5][5][0][0] = L3 + L4 - ReL5;
    Curvature_Higgs_L4[5][5][0][1] = -ImL5;
    Curvature_Higgs_L4[5][5][1][0] = -ImL5;
    Curvature_Higgs_L4[5][5][1][1] = L3 + L4 + ReL5;
    Curvature_Higgs_L4[5][5][2][2] = L3;
    Curvature_Higgs_L4[5][5][3][3] = L3;
    Curvature_Higgs_L4[5][5][4][4] = L2;
    Curvature_Higgs_L4[5][5][5][5] = 3 * L2;
    Curvature_Higgs_L4[5][5][6][6] = L2;
    Curvature_Higgs_L4[5][5][7][7] = L2;
    Curvature_Higgs_L4[5][5][8][8] = L8 / 0.2e1;
    Curvature_Higgs_L4[5][5][9][9] = L8 / 0.2e1;
    Curvature_Higgs_L4[5][6][0][2] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[5][6][0][3] = -L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[5][6][1][2] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[5][6][1][3] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[5][6][2][0] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[5][6][2][1] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[5][6][3][0] = -L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[5][6][3][1] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[5][6][5][6] = L2;
    Curvature_Higgs_L4[5][6][6][5] = L2;
    Curvature_Higgs_L4[5][7][0][2] = L4 / 0.2e1 - ReL5 / 0.2e1;
    Curvature_Higgs_L4[5][7][0][3] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[5][7][1][2] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[5][7][1][3] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[5][7][2][0] = L4 / 0.2e1 - ReL5 / 0.2e1;
    Curvature_Higgs_L4[5][7][2][1] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[5][7][3][0] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[5][7][3][1] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[5][7][5][7] = L2;
    Curvature_Higgs_L4[5][7][7][5] = L2;
    Curvature_Higgs_L4[5][8][5][8] = L8 / 0.2e1;
    Curvature_Higgs_L4[5][8][8][5] = L8 / 0.2e1;
    Curvature_Higgs_L4[5][9][5][9] = L8 / 0.2e1;
    Curvature_Higgs_L4[5][9][9][5] = L8 / 0.2e1;
    Curvature_Higgs_L4[6][0][0][6] = L3;
    Curvature_Higgs_L4[6][0][2][4] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[6][0][2][5] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[6][0][3][4] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[6][0][3][5] = -L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[6][0][4][2] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[6][0][4][3] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[6][0][5][2] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[6][0][5][3] = -L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[6][0][6][0] = L3;
    Curvature_Higgs_L4[6][1][1][6] = L3;
    Curvature_Higgs_L4[6][1][2][4] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[6][1][2][5] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[6][1][3][4] = L4 / 0.2e1 - ReL5 / 0.2e1;
    Curvature_Higgs_L4[6][1][3][5] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[6][1][4][2] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[6][1][4][3] = L4 / 0.2e1 - ReL5 / 0.2e1;
    Curvature_Higgs_L4[6][1][5][2] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[6][1][5][3] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[6][1][6][1] = L3;
    Curvature_Higgs_L4[6][2][0][4] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[6][2][0][5] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[6][2][1][4] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[6][2][1][5] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[6][2][2][6] = L3 + L4 + ReL5;
    Curvature_Higgs_L4[6][2][2][7] = -ImL5;
    Curvature_Higgs_L4[6][2][3][6] = ImL5;
    Curvature_Higgs_L4[6][2][3][7] = ReL5;
    Curvature_Higgs_L4[6][2][4][0] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[6][2][4][1] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[6][2][5][0] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[6][2][5][1] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[6][2][6][2] = L3 + L4 + ReL5;
    Curvature_Higgs_L4[6][2][6][3] = ImL5;
    Curvature_Higgs_L4[6][2][7][2] = -ImL5;
    Curvature_Higgs_L4[6][2][7][3] = ReL5;
    Curvature_Higgs_L4[6][3][0][4] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[6][3][0][5] = -L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[6][3][1][4] = L4 / 0.2e1 - ReL5 / 0.2e1;
    Curvature_Higgs_L4[6][3][1][5] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[6][3][2][6] = ImL5;
    Curvature_Higgs_L4[6][3][2][7] = ReL5;
    Curvature_Higgs_L4[6][3][3][6] = L3 + L4 - ReL5;
    Curvature_Higgs_L4[6][3][3][7] = ImL5;
    Curvature_Higgs_L4[6][3][4][0] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[6][3][4][1] = L4 / 0.2e1 - ReL5 / 0.2e1;
    Curvature_Higgs_L4[6][3][5][0] = -L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[6][3][5][1] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[6][3][6][2] = ImL5;
    Curvature_Higgs_L4[6][3][6][3] = L3 + L4 - ReL5;
    Curvature_Higgs_L4[6][3][7][2] = ReL5;
    Curvature_Higgs_L4[6][3][7][3] = ImL5;
    Curvature_Higgs_L4[6][4][0][2] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[6][4][0][3] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[6][4][1][2] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[6][4][1][3] = L4 / 0.2e1 - ReL5 / 0.2e1;
    Curvature_Higgs_L4[6][4][2][0] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[6][4][2][1] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[6][4][3][0] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[6][4][3][1] = L4 / 0.2e1 - ReL5 / 0.2e1;
    Curvature_Higgs_L4[6][4][4][6] = L2;
    Curvature_Higgs_L4[6][4][6][4] = L2;
    Curvature_Higgs_L4[6][5][0][2] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[6][5][0][3] = -L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[6][5][1][2] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[6][5][1][3] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[6][5][2][0] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[6][5][2][1] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[6][5][3][0] = -L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[6][5][3][1] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[6][5][5][6] = L2;
    Curvature_Higgs_L4[6][5][6][5] = L2;
    Curvature_Higgs_L4[6][6][0][0] = L3;
    Curvature_Higgs_L4[6][6][1][1] = L3;
    Curvature_Higgs_L4[6][6][2][2] = L3 + L4 + ReL5;
    Curvature_Higgs_L4[6][6][2][3] = ImL5;
    Curvature_Higgs_L4[6][6][3][2] = ImL5;
    Curvature_Higgs_L4[6][6][3][3] = L3 + L4 - ReL5;
    Curvature_Higgs_L4[6][6][4][4] = L2;
    Curvature_Higgs_L4[6][6][5][5] = L2;
    Curvature_Higgs_L4[6][6][6][6] = 3 * L2;
    Curvature_Higgs_L4[6][6][7][7] = L2;
    Curvature_Higgs_L4[6][6][8][8] = L8 / 0.2e1;
    Curvature_Higgs_L4[6][6][9][9] = L8 / 0.2e1;
    Curvature_Higgs_L4[6][7][2][2] = -ImL5;
    Curvature_Higgs_L4[6][7][2][3] = ReL5;
    Curvature_Higgs_L4[6][7][3][2] = ReL5;
    Curvature_Higgs_L4[6][7][3][3] = ImL5;
    Curvature_Higgs_L4[6][7][6][7] = L2;
    Curvature_Higgs_L4[6][7][7][6] = L2;
    Curvature_Higgs_L4[6][8][6][8] = L8 / 0.2e1;
    Curvature_Higgs_L4[6][8][8][6] = L8 / 0.2e1;
    Curvature_Higgs_L4[6][9][6][9] = L8 / 0.2e1;
    Curvature_Higgs_L4[6][9][9][6] = L8 / 0.2e1;
    Curvature_Higgs_L4[7][0][0][7] = L3;
    Curvature_Higgs_L4[7][0][2][4] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[7][0][2][5] = L4 / 0.2e1 - ReL5 / 0.2e1;
    Curvature_Higgs_L4[7][0][3][4] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[7][0][3][5] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[7][0][4][2] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[7][0][4][3] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[7][0][5][2] = L4 / 0.2e1 - ReL5 / 0.2e1;
    Curvature_Higgs_L4[7][0][5][3] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[7][0][7][0] = L3;
    Curvature_Higgs_L4[7][1][1][7] = L3;
    Curvature_Higgs_L4[7][1][2][4] = -L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[7][1][2][5] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[7][1][3][4] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[7][1][3][5] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[7][1][4][2] = -L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[7][1][4][3] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[7][1][5][2] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[7][1][5][3] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[7][1][7][1] = L3;
    Curvature_Higgs_L4[7][2][0][4] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[7][2][0][5] = L4 / 0.2e1 - ReL5 / 0.2e1;
    Curvature_Higgs_L4[7][2][1][4] = -L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[7][2][1][5] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[7][2][2][6] = -ImL5;
    Curvature_Higgs_L4[7][2][2][7] = L3 + L4 - ReL5;
    Curvature_Higgs_L4[7][2][3][6] = ReL5;
    Curvature_Higgs_L4[7][2][3][7] = -ImL5;
    Curvature_Higgs_L4[7][2][4][0] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[7][2][4][1] = -L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[7][2][5][0] = L4 / 0.2e1 - ReL5 / 0.2e1;
    Curvature_Higgs_L4[7][2][5][1] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[7][2][6][2] = -ImL5;
    Curvature_Higgs_L4[7][2][6][3] = ReL5;
    Curvature_Higgs_L4[7][2][7][2] = L3 + L4 - ReL5;
    Curvature_Higgs_L4[7][2][7][3] = -ImL5;
    Curvature_Higgs_L4[7][3][0][4] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[7][3][0][5] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[7][3][1][4] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[7][3][1][5] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[7][3][2][6] = ReL5;
    Curvature_Higgs_L4[7][3][2][7] = -ImL5;
    Curvature_Higgs_L4[7][3][3][6] = ImL5;
    Curvature_Higgs_L4[7][3][3][7] = L3 + L4 + ReL5;
    Curvature_Higgs_L4[7][3][4][0] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[7][3][4][1] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[7][3][5][0] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[7][3][5][1] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[7][3][6][2] = ReL5;
    Curvature_Higgs_L4[7][3][6][3] = ImL5;
    Curvature_Higgs_L4[7][3][7][2] = -ImL5;
    Curvature_Higgs_L4[7][3][7][3] = L3 + L4 + ReL5;
    Curvature_Higgs_L4[7][4][0][2] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[7][4][0][3] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[7][4][1][2] = -L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[7][4][1][3] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[7][4][2][0] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[7][4][2][1] = -L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[7][4][3][0] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[7][4][3][1] = ImL5 / 0.2e1;
    Curvature_Higgs_L4[7][4][4][7] = L2;
    Curvature_Higgs_L4[7][4][7][4] = L2;
    Curvature_Higgs_L4[7][5][0][2] = L4 / 0.2e1 - ReL5 / 0.2e1;
    Curvature_Higgs_L4[7][5][0][3] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[7][5][1][2] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[7][5][1][3] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[7][5][2][0] = L4 / 0.2e1 - ReL5 / 0.2e1;
    Curvature_Higgs_L4[7][5][2][1] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[7][5][3][0] = -ImL5 / 0.2e1;
    Curvature_Higgs_L4[7][5][3][1] = L4 / 0.2e1 + ReL5 / 0.2e1;
    Curvature_Higgs_L4[7][5][5][7] = L2;
    Curvature_Higgs_L4[7][5][7][5] = L2;
    Curvature_Higgs_L4[7][6][2][2] = -ImL5;
    Curvature_Higgs_L4[7][6][2][3] = ReL5;
    Curvature_Higgs_L4[7][6][3][2] = ReL5;
    Curvature_Higgs_L4[7][6][3][3] = ImL5;
    Curvature_Higgs_L4[7][6][6][7] = L2;
    Curvature_Higgs_L4[7][6][7][6] = L2;
    Curvature_Higgs_L4[7][7][0][0] = L3;
    Curvature_Higgs_L4[7][7][1][1] = L3;
    Curvature_Higgs_L4[7][7][2][2] = L3 + L4 - ReL5;
    Curvature_Higgs_L4[7][7][2][3] = -ImL5;
    Curvature_Higgs_L4[7][7][3][2] = -ImL5;
    Curvature_Higgs_L4[7][7][3][3] = L3 + L4 + ReL5;
    Curvature_Higgs_L4[7][7][4][4] = L2;
    Curvature_Higgs_L4[7][7][5][5] = L2;
    Curvature_Higgs_L4[7][7][6][6] = L2;
    Curvature_Higgs_L4[7][7][7][7] = 3 * L2;
    Curvature_Higgs_L4[7][7][8][8] = L8 / 0.2e1;
    Curvature_Higgs_L4[7][7][9][9] = L8 / 0.2e1;
    Curvature_Higgs_L4[7][8][7][8] = L8 / 0.2e1;
    Curvature_Higgs_L4[7][8][8][7] = L8 / 0.2e1;
    Curvature_Higgs_L4[7][9][7][9] = L8 / 0.2e1;
    Curvature_Higgs_L4[7][9][9][7] = L8 / 0.2e1;
    Curvature_Higgs_L4[8][0][0][8] = L7 / 0.2e1;
    Curvature_Higgs_L4[8][0][8][0] = L7 / 0.2e1;
    Curvature_Higgs_L4[8][1][1][8] = L7 / 0.2e1;
    Curvature_Higgs_L4[8][1][8][1] = L7 / 0.2e1;
    Curvature_Higgs_L4[8][2][2][8] = L7 / 0.2e1;
    Curvature_Higgs_L4[8][2][8][2] = L7 / 0.2e1;
    Curvature_Higgs_L4[8][3][3][8] = L7 / 0.2e1;
    Curvature_Higgs_L4[8][3][8][3] = L7 / 0.2e1;
    Curvature_Higgs_L4[8][4][4][8] = L8 / 0.2e1;
    Curvature_Higgs_L4[8][4][8][4] = L8 / 0.2e1;
    Curvature_Higgs_L4[8][5][5][8] = L8 / 0.2e1;
    Curvature_Higgs_L4[8][5][8][5] = L8 / 0.2e1;
    Curvature_Higgs_L4[8][6][6][8] = L8 / 0.2e1;
    Curvature_Higgs_L4[8][6][8][6] = L8 / 0.2e1;
    Curvature_Higgs_L4[8][7][7][8] = L8 / 0.2e1;
    Curvature_Higgs_L4[8][7][8][7] = L8 / 0.2e1;
    Curvature_Higgs_L4[8][8][0][0] = L7 / 0.2e1;
    Curvature_Higgs_L4[8][8][1][1] = L7 / 0.2e1;
    Curvature_Higgs_L4[8][8][2][2] = L7 / 0.2e1;
    Curvature_Higgs_L4[8][8][3][3] = L7 / 0.2e1;
    Curvature_Higgs_L4[8][8][4][4] = L8 / 0.2e1;
    Curvature_Higgs_L4[8][8][5][5] = L8 / 0.2e1;
    Curvature_Higgs_L4[8][8][6][6] = L8 / 0.2e1;
    Curvature_Higgs_L4[8][8][7][7] = L8 / 0.2e1;
    Curvature_Higgs_L4[8][8][8][8] = 0.3e1 / 0.4e1 * L6;
    Curvature_Higgs_L4[8][8][9][9] = L6 / 0.4e1;
    Curvature_Higgs_L4[8][9][8][9] = L6 / 0.4e1;
    Curvature_Higgs_L4[8][9][9][8] = L6 / 0.4e1;
    Curvature_Higgs_L4[9][0][0][9] = L7 / 0.2e1;
    Curvature_Higgs_L4[9][0][9][0] = L7 / 0.2e1;
    Curvature_Higgs_L4[9][1][1][9] = L7 / 0.2e1;
    Curvature_Higgs_L4[9][1][9][1] = L7 / 0.2e1;
    Curvature_Higgs_L4[9][2][2][9] = L7 / 0.2e1;
    Curvature_Higgs_L4[9][2][9][2] = L7 / 0.2e1;
    Curvature_Higgs_L4[9][3][3][9] = L7 / 0.2e1;
    Curvature_Higgs_L4[9][3][9][3] = L7 / 0.2e1;
    Curvature_Higgs_L4[9][4][4][9] = L8 / 0.2e1;
    Curvature_Higgs_L4[9][4][9][4] = L8 / 0.2e1;
    Curvature_Higgs_L4[9][5][5][9] = L8 / 0.2e1;
    Curvature_Higgs_L4[9][5][9][5] = L8 / 0.2e1;
    Curvature_Higgs_L4[9][6][6][9] = L8 / 0.2e1;
    Curvature_Higgs_L4[9][6][9][6] = L8 / 0.2e1;
    Curvature_Higgs_L4[9][7][7][9] = L8 / 0.2e1;
    Curvature_Higgs_L4[9][7][9][7] = L8 / 0.2e1;
    Curvature_Higgs_L4[9][8][8][9] = L6 / 0.4e1;
    Curvature_Higgs_L4[9][8][9][8] = L6 / 0.4e1;
    Curvature_Higgs_L4[9][9][0][0] = L7 / 0.2e1;
    Curvature_Higgs_L4[9][9][1][1] = L7 / 0.2e1;
    Curvature_Higgs_L4[9][9][2][2] = L7 / 0.2e1;
    Curvature_Higgs_L4[9][9][3][3] = L7 / 0.2e1;
    Curvature_Higgs_L4[9][9][4][4] = L8 / 0.2e1;
    Curvature_Higgs_L4[9][9][5][5] = L8 / 0.2e1;
    Curvature_Higgs_L4[9][9][6][6] = L8 / 0.2e1;
    Curvature_Higgs_L4[9][9][7][7] = L8 / 0.2e1;
    Curvature_Higgs_L4[9][9][8][8] = L6 / 0.4e1;
    Curvature_Higgs_L4[9][9][9][9] = 0.3e1 / 0.4e1 * L6;


    //    Curvature_Gauge_G2H2

    Curvature_Gauge_G2H2[0][0][0][0] = C_g * C_g / 0.2e1;
    Curvature_Gauge_G2H2[0][0][1][1] = C_g * C_g / 0.2e1;
    Curvature_Gauge_G2H2[0][0][2][2] = C_g * C_g / 0.2e1;
    Curvature_Gauge_G2H2[0][0][3][3] = C_g * C_g / 0.2e1;
    Curvature_Gauge_G2H2[0][0][4][4] = C_g * C_g / 0.2e1;
    Curvature_Gauge_G2H2[0][0][5][5] = C_g * C_g / 0.2e1;
    Curvature_Gauge_G2H2[0][0][6][6] = C_g * C_g / 0.2e1;
    Curvature_Gauge_G2H2[0][0][7][7] = C_g * C_g / 0.2e1;
    Curvature_Gauge_G2H2[0][3][0][2] = C_gs * C_g / 0.2e1;
    Curvature_Gauge_G2H2[0][3][1][3] = C_gs * C_g / 0.2e1;
    Curvature_Gauge_G2H2[0][3][2][0] = C_gs * C_g / 0.2e1;
    Curvature_Gauge_G2H2[0][3][3][1] = C_gs * C_g / 0.2e1;
    Curvature_Gauge_G2H2[0][3][4][6] = C_gs * C_g / 0.2e1;
    Curvature_Gauge_G2H2[0][3][5][7] = C_gs * C_g / 0.2e1;
    Curvature_Gauge_G2H2[0][3][6][4] = C_gs * C_g / 0.2e1;
    Curvature_Gauge_G2H2[0][3][7][5] = C_gs * C_g / 0.2e1;
    Curvature_Gauge_G2H2[1][1][0][0] = C_g * C_g / 0.2e1;
    Curvature_Gauge_G2H2[1][1][1][1] = C_g * C_g / 0.2e1;
    Curvature_Gauge_G2H2[1][1][2][2] = C_g * C_g / 0.2e1;
    Curvature_Gauge_G2H2[1][1][3][3] = C_g * C_g / 0.2e1;
    Curvature_Gauge_G2H2[1][1][4][4] = C_g * C_g / 0.2e1;
    Curvature_Gauge_G2H2[1][1][5][5] = C_g * C_g / 0.2e1;
    Curvature_Gauge_G2H2[1][1][6][6] = C_g * C_g / 0.2e1;
    Curvature_Gauge_G2H2[1][1][7][7] = C_g * C_g / 0.2e1;
    Curvature_Gauge_G2H2[1][3][0][3] = C_gs * C_g / 0.2e1;
    Curvature_Gauge_G2H2[1][3][1][2] = -C_gs * C_g / 0.2e1;
    Curvature_Gauge_G2H2[1][3][2][1] = -C_gs * C_g / 0.2e1;
    Curvature_Gauge_G2H2[1][3][3][0] = C_gs * C_g / 0.2e1;
    Curvature_Gauge_G2H2[1][3][4][7] = C_gs * C_g / 0.2e1;
    Curvature_Gauge_G2H2[1][3][5][6] = -C_gs * C_g / 0.2e1;
    Curvature_Gauge_G2H2[1][3][6][5] = -C_gs * C_g / 0.2e1;
    Curvature_Gauge_G2H2[1][3][7][4] = C_gs * C_g / 0.2e1;
    Curvature_Gauge_G2H2[2][2][0][0] = C_g * C_g / 0.2e1;
    Curvature_Gauge_G2H2[2][2][1][1] = C_g * C_g / 0.2e1;
    Curvature_Gauge_G2H2[2][2][2][2] = C_g * C_g / 0.2e1;
    Curvature_Gauge_G2H2[2][2][3][3] = C_g * C_g / 0.2e1;
    Curvature_Gauge_G2H2[2][2][4][4] = C_g * C_g / 0.2e1;
    Curvature_Gauge_G2H2[2][2][5][5] = C_g * C_g / 0.2e1;
    Curvature_Gauge_G2H2[2][2][6][6] = C_g * C_g / 0.2e1;
    Curvature_Gauge_G2H2[2][2][7][7] = C_g * C_g / 0.2e1;
    Curvature_Gauge_G2H2[2][3][0][0] = C_gs * C_g / 0.2e1;
    Curvature_Gauge_G2H2[2][3][1][1] = C_gs * C_g / 0.2e1;
    Curvature_Gauge_G2H2[2][3][2][2] = -C_gs * C_g / 0.2e1;
    Curvature_Gauge_G2H2[2][3][3][3] = -C_gs * C_g / 0.2e1;
    Curvature_Gauge_G2H2[2][3][4][4] = C_gs * C_g / 0.2e1;
    Curvature_Gauge_G2H2[2][3][5][5] = C_gs * C_g / 0.2e1;
    Curvature_Gauge_G2H2[2][3][6][6] = -C_gs * C_g / 0.2e1;
    Curvature_Gauge_G2H2[2][3][7][7] = -C_gs * C_g / 0.2e1;
    Curvature_Gauge_G2H2[3][0][0][2] = C_gs * C_g / 0.2e1;
    Curvature_Gauge_G2H2[3][0][1][3] = C_gs * C_g / 0.2e1;
    Curvature_Gauge_G2H2[3][0][2][0] = C_gs * C_g / 0.2e1;
    Curvature_Gauge_G2H2[3][0][3][1] = C_gs * C_g / 0.2e1;
    Curvature_Gauge_G2H2[3][0][4][6] = C_gs * C_g / 0.2e1;
    Curvature_Gauge_G2H2[3][0][5][7] = C_gs * C_g / 0.2e1;
    Curvature_Gauge_G2H2[3][0][6][4] = C_gs * C_g / 0.2e1;
    Curvature_Gauge_G2H2[3][0][7][5] = C_gs * C_g / 0.2e1;
    Curvature_Gauge_G2H2[3][1][0][3] = C_gs * C_g / 0.2e1;
    Curvature_Gauge_G2H2[3][1][1][2] = -C_gs * C_g / 0.2e1;
    Curvature_Gauge_G2H2[3][1][2][1] = -C_gs * C_g / 0.2e1;
    Curvature_Gauge_G2H2[3][1][3][0] = C_gs * C_g / 0.2e1;
    Curvature_Gauge_G2H2[3][1][4][7] = C_gs * C_g / 0.2e1;
    Curvature_Gauge_G2H2[3][1][5][6] = -C_gs * C_g / 0.2e1;
    Curvature_Gauge_G2H2[3][1][6][5] = -C_gs * C_g / 0.2e1;
    Curvature_Gauge_G2H2[3][1][7][4] = C_gs * C_g / 0.2e1;
    Curvature_Gauge_G2H2[3][2][0][0] = C_gs * C_g / 0.2e1;
    Curvature_Gauge_G2H2[3][2][1][1] = C_gs * C_g / 0.2e1;
    Curvature_Gauge_G2H2[3][2][2][2] = -C_gs * C_g / 0.2e1;
    Curvature_Gauge_G2H2[3][2][3][3] = -C_gs * C_g / 0.2e1;
    Curvature_Gauge_G2H2[3][2][4][4] = C_gs * C_g / 0.2e1;
    Curvature_Gauge_G2H2[3][2][5][5] = C_gs * C_g / 0.2e1;
    Curvature_Gauge_G2H2[3][2][6][6] = -C_gs * C_g / 0.2e1;
    Curvature_Gauge_G2H2[3][2][7][7] = -C_gs * C_g / 0.2e1;
    Curvature_Gauge_G2H2[3][3][0][0] = C_gs * C_gs / 0.2e1;
    Curvature_Gauge_G2H2[3][3][1][1] = C_gs * C_gs / 0.2e1;
    Curvature_Gauge_G2H2[3][3][2][2] = C_gs * C_gs / 0.2e1;
    Curvature_Gauge_G2H2[3][3][3][3] = C_gs * C_gs / 0.2e1;
    Curvature_Gauge_G2H2[3][3][4][4] = C_gs * C_gs / 0.2e1;
    Curvature_Gauge_G2H2[3][3][5][5] = C_gs * C_gs / 0.2e1;
    Curvature_Gauge_G2H2[3][3][6][6] = C_gs * C_gs / 0.2e1;
    Curvature_Gauge_G2H2[3][3][7][7] = C_gs * C_gs / 0.2e1;

    std::complex<double> II(0,1);

    // different types of 2HDM
    double vL = v2;
    double vd = v2;
    if(Type == 2) { vL = v1; vd = v1;}
    else if(Type == 3) vL = v1;
    else if(Type == 4) vd = v1;

    // set up coupling matrices

    MatrixXcd YIJ_zeta2(NQuarks,NQuarks), YIJ_eta2(NQuarks,NQuarks), YIJ_rho2(NQuarks,NQuarks), YIJ_psi2(NQuarks,NQuarks);

    MatrixXcd YIJ_zetaDown(NQuarks,NQuarks), YIJ_etaDown(NQuarks,NQuarks), YIJ_rhoDown(NQuarks,NQuarks), YIJ_psiDown(NQuarks,NQuarks);

    MatrixXcd YIJ_zetaLep(NLepton,NLepton), YIJ_etaLep(NLepton,NLepton), YIJ_rhoLep(NLepton,NLepton), YIJ_psiLep(NLepton,NLepton);

    YIJ_zeta2 = MatrixXcd::Zero(NQuarks,NQuarks);
    YIJ_eta2 = MatrixXcd::Zero(NQuarks,NQuarks);
    YIJ_rho2 = MatrixXcd::Zero(NQuarks,NQuarks);
    YIJ_psi2 = MatrixXcd::Zero(NQuarks,NQuarks);

    YIJ_zetaDown = MatrixXcd::Zero(NQuarks,NQuarks);
    YIJ_etaDown = MatrixXcd::Zero(NQuarks,NQuarks);
    YIJ_rhoDown = MatrixXcd::Zero(NQuarks,NQuarks);
    YIJ_psiDown = MatrixXcd::Zero(NQuarks,NQuarks);

    YIJ_zetaLep = MatrixXcd::Zero(NLepton,NLepton);
    YIJ_etaLep = MatrixXcd::Zero(NLepton,NLepton);
    YIJ_rhoLep = MatrixXcd::Zero(NLepton,NLepton);
    YIJ_psiLep = MatrixXcd::Zero(NLepton,NLepton);

    std::complex<double> V11,V12,V13,V21,V22,V23,V31,V32,V33;
    V11 = C_Vud;
    V12 = C_Vus;
    V13 = C_Vub;
    V21 = C_Vcd;
    V22 = C_Vcs;
    V23 = C_Vcb;
    V31 = C_Vtd;
    V32 = C_Vts;
    V33 = C_Vtb;

    //leptons
    YIJ_zetaLep(1,6) = 0.1e1 / vL * C_MassElectron;
    YIJ_zetaLep(3,7) = 0.1e1 / vL * C_MassMu;
    YIJ_zetaLep(5,8) = 0.1e1 / vL * C_MassTau;
    YIJ_zetaLep(6,1) = 0.1e1 / vL * C_MassElectron;
    YIJ_zetaLep(7,3) = 0.1e1 / vL * C_MassMu;
    YIJ_zetaLep(8,5) = 0.1e1 / vL * C_MassTau;
    YIJ_etaLep(1,6) = II / vL * C_MassElectron;
    YIJ_etaLep(3,7) = II / vL * C_MassMu;
    YIJ_etaLep(5,8) = II / vL * C_MassTau;
    YIJ_etaLep(6,1) = II / vL * C_MassElectron;
    YIJ_etaLep(7,3) = II / vL * C_MassMu;
    YIJ_etaLep(8,5) = II / vL * C_MassTau;
    YIJ_rhoLep(0,1) = 0.1e1 / vL * C_MassElectron;
    YIJ_rhoLep(1,0) = 0.1e1 / vL * C_MassElectron;
    YIJ_rhoLep(2,3) = 0.1e1 / vL * C_MassMu;
    YIJ_rhoLep(3,2) = 0.1e1 / vL * C_MassMu;
    YIJ_rhoLep(4,5) = 0.1e1 / vL * C_MassTau;
    YIJ_rhoLep(5,4) = 0.1e1 / vL * C_MassTau;
    YIJ_psiLep(0,1) = II / vL * C_MassElectron;
    YIJ_psiLep(1,0) = II / vL * C_MassElectron;
    YIJ_psiLep(2,3) = II / vL * C_MassMu;
    YIJ_psiLep(3,2) = II / vL * C_MassMu;
    YIJ_psiLep(4,5) = II / vL * C_MassTau;
    YIJ_psiLep(5,4) = II / vL * C_MassTau;

    YIJ_zetaDown(3,6) = V11 / vd * C_MassDown;
    YIJ_zetaDown(3,7) = V21 / vd * C_MassDown;
    YIJ_zetaDown(3,8) = V31 / vd * C_MassDown;
    YIJ_zetaDown(4,6) = V12 / vd * C_MassStrange;
    YIJ_zetaDown(4,7) = V22 / vd * C_MassStrange;
    YIJ_zetaDown(4,8) = V32 / vd * C_MassStrange;
    YIJ_zetaDown(5,6) = V13 / vd * C_MassBottom;
    YIJ_zetaDown(5,7) = V23 / vd * C_MassBottom;
    YIJ_zetaDown(5,8) = V33 / vd * C_MassBottom;
    YIJ_zetaDown(6,3) = V11 / vd * C_MassDown;
    YIJ_zetaDown(6,4) = V12 / vd * C_MassStrange;
    YIJ_zetaDown(6,5) = V13 / vd * C_MassBottom;
    YIJ_zetaDown(7,3) = V21 / vd * C_MassDown;
    YIJ_zetaDown(7,4) = V22 / vd * C_MassStrange;
    YIJ_zetaDown(7,5) = V23 / vd * C_MassBottom;
    YIJ_zetaDown(8,3) = V31 / vd * C_MassDown;
    YIJ_zetaDown(8,4) = V32 / vd * C_MassStrange;
    YIJ_zetaDown(8,5) = V33 / vd * C_MassBottom;
    YIJ_etaDown(3,6) = II * V11 / vd * C_MassDown;
    YIJ_etaDown(3,7) = II * V21 / vd * C_MassDown;
    YIJ_etaDown(3,8) = II * V31 / vd * C_MassDown;
    YIJ_etaDown(4,6) = II * V12 / vd * C_MassStrange;
    YIJ_etaDown(4,7) = II * V22 / vd * C_MassStrange;
    YIJ_etaDown(4,8) = II * V32 / vd * C_MassStrange;
    YIJ_etaDown(5,6) = II * V13 / vd * C_MassBottom;
    YIJ_etaDown(5,7) = II * V23 / vd * C_MassBottom;
    YIJ_etaDown(5,8) = II * V33 / vd * C_MassBottom;
    YIJ_etaDown(6,3) = II * V11 / vd * C_MassDown;
    YIJ_etaDown(6,4) = II * V12 / vd * C_MassStrange;
    YIJ_etaDown(6,5) = II * V13 / vd * C_MassBottom;
    YIJ_etaDown(7,3) = II * V21 / vd * C_MassDown;
    YIJ_etaDown(7,4) = II * V22 / vd * C_MassStrange;
    YIJ_etaDown(7,5) = II * V23 / vd * C_MassBottom;
    YIJ_etaDown(8,3) = II * V31 / vd * C_MassDown;
    YIJ_etaDown(8,4) = II * V32 / vd * C_MassStrange;
    YIJ_etaDown(8,5) = II * V33 / vd * C_MassBottom;
    YIJ_rhoDown(3,9) = 0.1e1 / vd * C_MassDown;
    YIJ_rhoDown(4,10) = 0.1e1 / vd * C_MassStrange;
    YIJ_rhoDown(5,11) = 0.1e1 / vd * C_MassBottom;
    YIJ_rhoDown(9,3) = 0.1e1 / vd * C_MassDown;
    YIJ_rhoDown(10,4) = 0.1e1 / vd * C_MassStrange;
    YIJ_rhoDown(11,5) = 0.1e1 / vd * C_MassBottom;
    YIJ_psiDown(3,9) = II / vd * C_MassDown;
    YIJ_psiDown(4,10) = II / vd * C_MassStrange;
    YIJ_psiDown(5,11) = II / vd * C_MassBottom;
    YIJ_psiDown(9,3) = II / vd * C_MassDown;
    YIJ_psiDown(10,4) = II / vd * C_MassStrange;
    YIJ_psiDown(11,5) = II / vd * C_MassBottom;

    //REDO THIS:
    YIJ_zeta2(0,9) = -0.1e1 / v2 * C_MassUp * conj(V11);
    YIJ_zeta2(0,10) = -0.1e1 / v2 * C_MassUp * conj(V12);
    YIJ_zeta2(0,11) = -0.1e1 / v2 * C_MassUp * conj(V13);
    YIJ_zeta2(1,9) = -0.1e1 / v2 * C_MassCharm * conj(V21);
    YIJ_zeta2(1,10) = -0.1e1 / v2 * C_MassCharm * conj(V22);
    YIJ_zeta2(1,11) = -0.1e1 / v2 * C_MassCharm * conj(V23);
    YIJ_zeta2(2,9) = -0.1e1 / v2 * C_MassTop * conj(V31);
    YIJ_zeta2(2,10) = -0.1e1 / v2 * C_MassTop * conj(V32);
    YIJ_zeta2(2,11) = -0.1e1 / v2 * C_MassTop * conj(V33);
    YIJ_zeta2(9,0) = -0.1e1 / v2 * C_MassUp * conj(V11);
    YIJ_zeta2(9,1) = -0.1e1 / v2 * C_MassCharm * conj(V21);
    YIJ_zeta2(9,2) = -0.1e1 / v2 * C_MassTop * conj(V31);
    YIJ_zeta2(10,0) = -0.1e1 / v2 * C_MassUp * conj(V12);
    YIJ_zeta2(10,1) = -0.1e1 / v2 * C_MassCharm * conj(V22);
    YIJ_zeta2(10,2) = -0.1e1 / v2 * C_MassTop * conj(V32);
    YIJ_zeta2(11,0) = -0.1e1 / v2 * C_MassUp * conj(V13);
    YIJ_zeta2(11,1) = -0.1e1 / v2 * C_MassCharm * conj(V23);
    YIJ_zeta2(11,2) = -0.1e1 / v2 * C_MassTop * conj(V33);
    YIJ_eta2(0,9) = II / v2 * C_MassUp * conj(V11);
    YIJ_eta2(0,10) = II / v2 * C_MassUp * conj(V12);
    YIJ_eta2(0,11) = II / v2 * C_MassUp * conj(V13);
    YIJ_eta2(1,9) = II / v2 * C_MassCharm * conj(V21);
    YIJ_eta2(1,10) = II / v2 * C_MassCharm * conj(V22);
    YIJ_eta2(1,11) = II / v2 * C_MassCharm * conj(V23);
    YIJ_eta2(2,9) = II / v2 * C_MassTop * conj(V31);
    YIJ_eta2(2,10) = II / v2 * C_MassTop * conj(V32);
    YIJ_eta2(2,11) = II / v2 * C_MassTop * conj(V33);
    YIJ_eta2(9,0) = II / v2 * C_MassUp * conj(V11);
    YIJ_eta2(9,1) = II / v2 * C_MassCharm * conj(V21);
    YIJ_eta2(9,2) = II / v2 * C_MassTop * conj(V31);
    YIJ_eta2(10,0) = II / v2 * C_MassUp * conj(V12);
    YIJ_eta2(10,1) = II / v2 * C_MassCharm * conj(V22);
    YIJ_eta2(10,2) = II / v2 * C_MassTop * conj(V32);
    YIJ_eta2(11,0) = II / v2 * C_MassUp * conj(V13);
    YIJ_eta2(11,1) = II / v2 * C_MassCharm * conj(V23);
    YIJ_eta2(11,2) = II / v2 * C_MassTop * conj(V33);
    YIJ_rho2(0,6) = 0.1e1 / v2 * C_MassUp;
    YIJ_rho2(1,7) = 0.1e1 / v2 * C_MassCharm;
    YIJ_rho2(2,8) = 0.1e1 / v2 * C_MassTop;
    YIJ_rho2(6,0) = 0.1e1 / v2 * C_MassUp;
    YIJ_rho2(7,1) = 0.1e1 / v2 * C_MassCharm;
    YIJ_rho2(8,2) = 0.1e1 / v2 * C_MassTop;
    YIJ_psi2(0,6) = -II / v2 * C_MassUp;
    YIJ_psi2(1,7) = -II / v2 * C_MassCharm;
    YIJ_psi2(2,8) = -II / v2 * C_MassTop;
    YIJ_psi2(6,0) = -II / v2 * C_MassUp;
    YIJ_psi2(7,1) = -II / v2 * C_MassCharm;
    YIJ_psi2(8,2) = -II / v2 * C_MassTop;








////// CHECK THIS!

    for(size_t i=0;i<NLepton;i++)
	{
		for(size_t j=0;j<NLepton;j++)
		{
			if(Type == 1 or Type == 4)
			{
                Curvature_Lepton_F2H1[i][j][4] = YIJ_zetaLep(i,j);
                Curvature_Lepton_F2H1[i][j][5] = YIJ_etaLep(i,j);
                Curvature_Lepton_F2H1[i][j][6] = YIJ_rhoLep(i,j);
                Curvature_Lepton_F2H1[i][j][7] = YIJ_psiLep(i,j);
			}
            //type 2 and type 3
			else{
                Curvature_Lepton_F2H1[i][j][0] = YIJ_zetaLep(i,j);
                Curvature_Lepton_F2H1[i][j][1] = YIJ_etaLep(i,j);
                Curvature_Lepton_F2H1[i][j][2] = YIJ_rhoLep(i,j);
                Curvature_Lepton_F2H1[i][j][3] = YIJ_psiLep(i,j);
			}
		}
	}

	for(size_t i=0;i<NQuarks;i++)
	{
		for(size_t j=0;j<NQuarks;j++)
		{
            Curvature_Quark_F2H1[i][j][4] = YIJ_zeta2(i,j);
            Curvature_Quark_F2H1[i][j][5] = YIJ_eta2(i,j);
            Curvature_Quark_F2H1[i][j][6] = YIJ_rho2(i,j);
            Curvature_Quark_F2H1[i][j][7] = YIJ_psi2(i,j);

			if(Type == 1 or Type == 3)
			{
                Curvature_Quark_F2H1[i][j][4] += YIJ_zetaDown(i,j);
                Curvature_Quark_F2H1[i][j][5] += YIJ_etaDown(i,j);
                Curvature_Quark_F2H1[i][j][6] += YIJ_rhoDown(i,j);
				Curvature_Quark_F2H1[i][j][7] += YIJ_psiDown(i,j);
			}
			else{
				Curvature_Quark_F2H1[i][j][0] += YIJ_zetaDown(i,j);
                Curvature_Quark_F2H1[i][j][1] += YIJ_etaDown(i,j);
                Curvature_Quark_F2H1[i][j][2] += YIJ_rhoDown(i,j);
                Curvature_Quark_F2H1[i][j][3] += YIJ_psiDown(i,j);
			}
		}
	}



}






bool Class_Potential_CN2HDM::CalculateDebyeSimplified(){
  return false;
  //optional
  // set if using extending fermion or gauge sectors - generic formula only works if just an extended Higgs sector.
  // Plan for code for this to work even if you have extended gauge and fermion sectors - just to work from tree level lagrangian.
  /*
   * Use this function if you calculated the Debye corrections to the Higgs mass matrix and implement
   * your formula here and return true. The vector is given by DebyeHiggs[NHiggs][NHiggs]
   */
}

bool Class_Potential_CN2HDM::CalculateDebyeGaugeSimplified()
{
    bool Debug = false;
  if(Debug) std::cout << "Debug turned on in Class_Potential_CN2HDM :: " << __func__ << std::endl;
//optional
  /*
     * Use this function if you calculated the Debye corrections to the gauge mass matrix and implement
     * your formula here and return true. The vector is given by DebyeGauge[NGauge][NGauge]
     */


  return false;
}
double Class_Potential_CN2HDM::VTreeSimplified(const std::vector<double>& v) const
{
    //optional
    (void )v;

	return 0;
}

double Class_Potential_CN2HDM::VCounterSimplified(const std::vector<double>& v) const
{
    (void) v;
	return 0;
}

void Class_Potential_CN2HDM::Debugging(const std::vector<double>& input, std::vector<double>& output) const
{
	(void) input;
	(void) output;

}

}
}
