/*
 * transport_equations.cpp
 *
 *  Copyright (C) 2020  Philipp Basler, Margarete Mühlleitner and Jonas Müller

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

#include <BSMPT/baryo_calculation/transport_equations.h>
#include <BSMPT/Kfactors/Kfactors.h>
#include <gsl/gsl_integration.h>
#include <BSMPT/minimizer/MinimizePlane.h>
#include <BSMPT/WallThickness/WallThicknessLib.h>
#include <BSMPT/models/ClassPotentialOrigin.h>
#include <BSMPT/Kfactors/KfactorsinterpolatedGSL.h>
#include <BSMPT/utility.h>

/**
 * @file
 */

namespace BSMPT{
namespace Baryo{

void GSL_integration_mubl::setSymmetricCPViolatingPhase(double Phase){
    symmetric_CP_violating_phase = Phase;
}

double GSL_integration_mubl::getSymmetricCPViolatingPhase() const{
    return symmetric_CP_violating_phase;
}

double GSL_integration_mubl::getSymmetricCPViolatingPhase_top() const{
    return TOP_symmetric_CP_violating_phase;
}

double GSL_integration_mubl::getSymmetricCPViolatingPhase_bot() const{
    return BOT_symmetric_CP_violating_phase;
}

double GSL_integration_mubl::getSymmetricCPViolatingPhase_tau() const{
    return TAU_symmetric_CP_violating_phase;
}


double GSL_integration_mubl::getBrokenCPViolatingPhase() const{
    return broken_CP_violating_phase;
}

double GSL_integration_mubl::getBrokenCPViolatingPhase_top() const{
    return TOP_broken_CP_violating_phase;
}
double GSL_integration_mubl::getBrokenCPViolatingPhase_bot() const{
    return BOT_broken_CP_violating_phase;
}
double GSL_integration_mubl::getBrokenCPViolatingPhase_tau() const{
    return TAU_broken_CP_violating_phase;
}

std::vector<double> GSL_integration_mubl::getVEVCritical() const{
    return vev_critical;
}
std::vector<double> GSL_integration_mubl::getVEVsym() const{
    return vev_symmetric;
}

void GSL_integration_mubl::setTC(double TC_in){
    TC = TC_in;
}

double GSL_integration_mubl::getTC() const{
    return TC;
}


void GSL_integration_mubl::setvw(double vw_in){
    vw = vw_in;
}

double GSL_integration_mubl::getvw() const{
    return vw;
}

void GSL_integration_mubl::setUseVelocityTransportEquations(bool in){
    UseVelocityTransportEquations = in;
}

bool GSL_integration_mubl::getUseVelocityTransportEquations() const{
    return UseVelocityTransportEquations;
}

std::shared_ptr<BSMPT::Class_Potential_Origin> GSL_integration_mubl::getModelPointer() const{
    return modelPointer;
}

double GSL_integration_mubl::getLW() const{
    return LW;
}

double GSL_integration_mubl::getZMAX() const{
    return zmax;
}

void GSL_integration_mubl::setpar(std::vector<double> inp){
    par = inp;
}
std::vector<double> GSL_integration_mubl::getpar(){
    return par;
}
TransportMethod GSL_integration_mubl::get_transport_method(){
    return transport_method;
}
void GSL_integration_mubl::set_transport_method(TransportMethod method){
    transport_method = method;
}
void GSL_integration_mubl::setZMAX(double zin,bool MultiplesOfLW=false){
    zmax = zin;
    if(MultiplesOfLW) zmax *= LW;
}
void GSL_integration_mubl::set_vev_sym_theta(std::vector<double> & vev_in){
    vev_sym_theta = vev_in;
}
std::vector<double> GSL_integration_mubl::get_vev_sym_theta() const{
    return vev_sym_theta;
}



void GSL_integration_mubl::init(const double& vw_input,
                                std::vector<double>& vev_critical_input,
                                std::vector<double>& vev_symmetric_input,
                                const double& TC_input,
                                std::shared_ptr<Class_Potential_Origin>& modelPointer_input,
                                const int& WhichMinimizer){

    vw = vw_input;
    TC = TC_input;
    vev_critical = vev_critical_input;
    vev_symmetric = vev_symmetric_input;
    modelPointer = modelPointer_input;
    setpar(modelPointer->get_parStored());
    UseVelocityTransportEquations = false;

    std::vector<double> vevCriticalNhiggs;
    vevCriticalNhiggs=modelPointer->MinimizeOrderVEV(vev_critical);
    vc = modelPointer->EWSBVEV(vevCriticalNhiggs);



    double CompareNorm = 0;
    for(auto x: modelPointer->get_vevTreeMin()) CompareNorm += std::pow(x,2);


    double TestNorm=0;
    for(std::size_t i =0 ;i<vev_critical.size();i++) TestNorm += std::pow(modelPointer->get_vevTreeMin(i) - vev_critical.at(i),2);
    if(TestNorm > CompareNorm)
    {
        for(auto symmetry : modelPointer->SignSymmetries){
            TestNorm = 0;
            for(std::size_t i =0 ;i<vev_critical.size();i++) {
                TestNorm += std::pow(modelPointer->get_vevTreeMin(i) - symmetry.at(i)*vev_critical.at(i),2);
            }
            if(TestNorm < CompareNorm){
                for(std::size_t i=0;i<vev_critical.size();i++) vev_critical.at(i) *= symmetry.at(i);
                break;
            }
        }
    }

    TestNorm=0;
    for(std::size_t i =0 ;i<vev_symmetric.size();i++) TestNorm += std::pow(modelPointer->get_vevTreeMin(i) - vev_symmetric.at(i),2);
    if(TestNorm > CompareNorm)
    {
        for(auto symmetry : modelPointer->SignSymmetries){
            TestNorm = 0;
            for(std::size_t i =0 ;i<vev_symmetric.size();i++) {
                TestNorm += std::pow(modelPointer->get_vevTreeMin(i) - symmetry.at(i)*vev_symmetric.at(i),2);
            }
            if(TestNorm < CompareNorm){
                for(std::size_t i=0;i<vev_symmetric.size();i++) vev_symmetric.at(i) *= symmetry.at(i);
                break;
            }
        }
    }

    LW = Wall::calculate_wall_thickness_plane(modelPointer,TC,vev_critical,vev_symmetric,WhichMinimizer);

    if(false)
    {
        double LW1D = Wall::calculate_wall_thickness_1D(modelPointer,TC,vev_critical,vev_symmetric);
        std::cout << "The 1D LW is given by " << LW1D << std::endl;
        std::cout << "The relative error to the plane LW is given by " << (LW-LW1D)/LW *100 << " %" << std::endl;
    }

    zmax = 4*LW;

    // Find minimum slightly before the symmetric phase to define the CP-violating phase in the symmetric minimum
    std::vector<double> basepoint;
    for(std::size_t i=0;i<vev_critical.size();i++)
    {
        basepoint.push_back(vev_symmetric.at(i) + 1e-2 * (vev_critical.at(i)-vev_symmetric.at(i)));
    }
    auto MinPlaneResult = Minimizer::MinimizePlane(basepoint,vev_symmetric,vev_critical,modelPointer,TC,WhichMinimizer);
    auto MinimumPlane = MinPlaneResult.Minimum;




    TestNorm = 0;
    for(std::size_t i =0 ;i<MinimumPlane.size();i++) TestNorm += std::pow(modelPointer->get_vevTreeMin(i) - MinimumPlane.at(i),2);
    if(TestNorm > CompareNorm)
    {
        for(auto symmetry : modelPointer->SignSymmetries){
            TestNorm = 0;
            for(std::size_t i =0 ;i<MinimumPlane.size();i++) {
                TestNorm += std::pow(modelPointer->get_vevTreeMin(i) - symmetry.at(i)*MinimumPlane.at(i),2);
            }
            if(TestNorm < CompareNorm){
                for(std::size_t i=0;i<MinimumPlane.size();i++) MinimumPlane.at(i) *= symmetry.at(i);
                break;
            }
        }
    }
    std::vector<double> TransformedVEV;
    std::vector<std::complex<double>> ComplexQuarkMassMinimumPlane;
    std::vector<std::complex<double>> ComplexLeptonMassMinimumPlane;
    TransformedVEV=modelPointer->MinimizeOrderVEV(MinimumPlane);
    ComplexQuarkMassMinimumPlane=modelPointer->QuarkMasses(TransformedVEV);
    ComplexLeptonMassMinimumPlane=modelPointer->LeptonMasses(TransformedVEV);
    symmetric_CP_violating_phase = std::arg(ComplexQuarkMassMinimumPlane.at(ComplexQuarkMassMinimumPlane.size()-1));
    TOP_symmetric_CP_violating_phase = std::arg(ComplexQuarkMassMinimumPlane.at(ComplexQuarkMassMinimumPlane.size()-1));
    BOT_symmetric_CP_violating_phase = std::arg(ComplexQuarkMassMinimumPlane.at(ComplexQuarkMassMinimumPlane.size()-3));
    TAU_symmetric_CP_violating_phase = std::arg(ComplexLeptonMassMinimumPlane.at(ComplexLeptonMassMinimumPlane.size()-1));
    if(std::abs(symmetric_CP_violating_phase) > 0.5*M_PI){
        symmetric_CP_violating_phase = std::arg(ComplexQuarkMassMinimumPlane.at(ComplexQuarkMassMinimumPlane.size()-2));
        TOP_symmetric_CP_violating_phase = std::arg(ComplexQuarkMassMinimumPlane.at(ComplexQuarkMassMinimumPlane.size()-2));
    }
    if(std::abs(BOT_symmetric_CP_violating_phase) > 0.5*M_PI){
        BOT_symmetric_CP_violating_phase = std::arg(ComplexQuarkMassMinimumPlane.at(ComplexQuarkMassMinimumPlane.size()-4));
    }
    if(std::abs(TAU_symmetric_CP_violating_phase) > 0.5*M_PI){
        TAU_symmetric_CP_violating_phase = std::arg(ComplexLeptonMassMinimumPlane.at(ComplexLeptonMassMinimumPlane.size()-2));
    }
    //Numerical Stability Check
    if(std::abs(TOP_symmetric_CP_violating_phase)<1e-14) TOP_symmetric_CP_violating_phase=0;
    if(std::abs(BOT_symmetric_CP_violating_phase)<1e-14) BOT_symmetric_CP_violating_phase=0;
    if(std::abs(TAU_symmetric_CP_violating_phase)<1e-14) TAU_symmetric_CP_violating_phase=0;


    std::vector<std::complex<double>> ComplexQuarkMassBroken;
    std::vector<std::complex<double>> ComplexLeptonMassBroken;
    TransformedVEV.clear();
    TransformedVEV=modelPointer->MinimizeOrderVEV(vev_critical);
    ComplexQuarkMassBroken=modelPointer->QuarkMasses(TransformedVEV);
    ComplexLeptonMassBroken=modelPointer->LeptonMasses(TransformedVEV);

    broken_CP_violating_phase = std::arg(ComplexQuarkMassBroken.at(ComplexQuarkMassBroken.size()-1));
    TOP_broken_CP_violating_phase = std::arg(ComplexQuarkMassBroken.at(ComplexQuarkMassBroken.size()-1));
    BOT_broken_CP_violating_phase = std::arg(ComplexQuarkMassBroken.at(ComplexQuarkMassBroken.size()-3));
    TAU_broken_CP_violating_phase = std::arg(ComplexLeptonMassBroken.at(ComplexLeptonMassBroken.size()-1));
    if(std::abs(broken_CP_violating_phase) > 0.5*M_PI){
        broken_CP_violating_phase = std::arg(ComplexQuarkMassBroken.at(ComplexQuarkMassBroken.size()-2));
        TOP_broken_CP_violating_phase = std::arg(ComplexQuarkMassBroken.at(ComplexQuarkMassBroken.size()-2));
    }
    if(std::abs(BOT_broken_CP_violating_phase) > 0.5*M_PI){
        BOT_broken_CP_violating_phase = std::arg(ComplexQuarkMassBroken.at(ComplexQuarkMassBroken.size()-4));
    }
    if(std::abs(TAU_broken_CP_violating_phase)> 0.5*M_PI){
        TAU_broken_CP_violating_phase = std::arg(ComplexLeptonMassBroken.at(ComplexLeptonMassBroken.size()-2));
    }
    //Numerical Stability Check
    if(std::abs(TOP_broken_CP_violating_phase)<1e-14) TOP_broken_CP_violating_phase=0;
    if(std::abs(BOT_broken_CP_violating_phase)<1e-14) BOT_broken_CP_violating_phase=0;
    if(std::abs(TAU_broken_CP_violating_phase)<1e-14) TAU_broken_CP_violating_phase=0;
    set_vev_sym_theta(MinimumPlane);
}

transport_equations::transport_equations(const struct GSL_integration_mubl& params)
    : UseVelocityTransportEquations{params.getUseVelocityTransportEquations()},
      vw{params.getvw()}, LW{params.getLW()},
      TC{params.getTC()},
      symmetric_CP_violating_phase{params.getSymmetricCPViolatingPhase()},
      broken_CP_violating_phase{params.getBrokenCPViolatingPhase()},
      modelPointer{params.getModelPointer()},
      vev_critical{modelPointer->MinimizeOrderVEV(params.getVEVCritical())}
{
    // TODO Auto-generated constructor stub

}

transport_equations::~transport_equations() {
    // TODO Auto-generated destructor stub
}





double transport_equations::get_W_mass(const std::vector<double>& vev, const double &T) const{
    std::vector<double> res;
    res=modelPointer->GaugeMassesSquared(vev,T);
    std::vector<double> nrepeat(modelPointer->get_NGauge());
    for(std::size_t i=0;i<modelPointer->get_NGauge();i++)
    {
        nrepeat[i]=0;
        for(std::size_t j=0;j<modelPointer->get_NGauge();j++){
            if(std::abs(res.at(i)- res.at(j)) <= 1e-5) nrepeat[i]++;
        }
    }

    for(int j=modelPointer->get_NGauge()-1;j>=0;j--){
        // std::cout << res.at(j) << std::endl;
        if(nrepeat[j]>1) {
            if(std::isnan(res.at(j))){
                std::string retmessage = "Nan found in ";
                retmessage+= __func__;
                throw std::runtime_error(retmessage);
            }
            return res.at(j);
        }
    }
    return 0;
}


std::vector<double> transport_equations::get_top_mass_and_derivative(const std::vector<double>& vev) const
{
    std::vector<double> restmp;
    std::vector<double> res;
    int postop=modelPointer->get_NQuarks() -1;
    int posdifftop=2*modelPointer->get_NQuarks() -1;
    double topmasssquared=-1;
    std::vector<double> topderivatives;
    for(std::size_t i=1;i<=modelPointer->get_NHiggs();i++)
    {
        restmp.clear();
        restmp=modelPointer->QuarkMassesSquared(vev,i);
        if(topmasssquared==-1) topmasssquared = restmp.at(postop);
        topderivatives.push_back(restmp.at(posdifftop));
    }
    res.push_back(topmasssquared);
    if(std::isnan(topmasssquared)){
        for(std::size_t i=0;i<vev.size();i++) std::cout << vev.at(i) << sep;
        std::cout << std::endl;
        std::string retmessage = "Nan found in ";
        retmessage+= __func__;
        throw std::runtime_error(retmessage);
    }
    for(std::size_t i=0;i<topderivatives.size();i++) {
        res.push_back(topderivatives.at(i));
        if(std::isnan(topderivatives.at(i))){
            std::string retmessage = "Nan found in ";
            retmessage+= __func__;
            retmessage+= " at deriv number ";
            retmessage+= std::to_string(i);
            for(std::size_t j=0;j<vev.size();j++) std::cout << vev.at(j) << sep;
            throw std::runtime_error(retmessage);
        }
    }
    return  res;
}


std::vector<double> transport_equations::calculate_vev(const double &z) const
{
    std::vector<double> vev;
    if(std::isnan(std::tanh(z/LW))){
        std::string retmessage = "Nan found in ";
        retmessage+= __func__;
        throw std::runtime_error(retmessage);
    }
    for(std::size_t i=0;i<modelPointer->get_NHiggs();i++){
        vev.push_back(0.5*vev_critical.at(i) *(1-std::tanh(z/LW)));
    }
    return vev;
}

std::vector<double> transport_equations::calculate_vev_derivative(const double &z) const
{
    std::vector<double> vev;
    double tanhv=std::tanh(z/LW);
    double diffval = -(1-tanhv*tanhv)/LW;
    for(std::size_t i=0;i<modelPointer->get_NHiggs();i++){
        vev.push_back(0.5*vev_critical.at(i) *diffval);
    }
    return vev;
}


double transport_equations::calculate_theta(const double& z, const int& diff) const
{
    double res = 0;
    double thetasym = symmetric_CP_violating_phase;
    if(modelPointer->get_Model() != ModelID::ModelIDs::C2HDM) std::cerr << "This is only programmed for the C2HDM" << std::endl;
    double  thetabrk = broken_CP_violating_phase;
    double difftheta = thetabrk - thetasym;
    double tanhv = std::tanh(z/LW);
    if(diff == 0){
        res = thetabrk - 0.5*difftheta*(1+tanhv);
    }
    else if(diff == 1){
        res = -difftheta/(2*LW) *(1-tanhv*tanhv);
    }
    else if(diff == 2){
        res = difftheta*tanhv/(LW*LW) *(1-tanhv*tanhv);
    }

    if(std::isnan(res)){
        std::string retmessage = "Nan found in ";
        retmessage+= __func__;
        throw std::runtime_error(retmessage);
    }

    double v3c=vev_critical.at(7), v2c=vev_critical.at(6), v1c=vev_critical.at(4);
    double tanbeta_temp_sq=(std::pow(v3c,2)+std::pow(v2c,2))/std::pow(v1c,2);

    if(UseTanBetaSuppression) res *= 1.0/(1+tanbeta_temp_sq); // Eq 27 of 0605242

    return res;
}


transport_equations &transport_equations::operator()(const state_type &x , state_type &dxdt , const double z ){
    std::cout << std::scientific;
    std::vector<double> topres;
    std::vector<double> vev,vevdiff;

    vev=calculate_vev(z);

    vevdiff=calculate_vev_derivative(z);

    topres=get_top_mass_and_derivative(vev);

    double mws = get_W_mass(vev,TC);

    double mtsquared = topres.at(0);

    double dmtsquared = 0;
    for(std::size_t i=0;i<vevdiff.size();i++){
        dmtsquared += vevdiff.at(i) * topres.at(i+1);
    }



    double K1t = Kfactors::K1fermion_normalized(mtsquared,TC);// K1_fermion_interp(mtsquared,TC,PyLoad)/norm1;
    double K2t = Kfactors::K2fermion_normalized(mtsquared,TC);// K2_fermion_interp(mtsquared,TC,PyLoad)/norm1;
    double K4t = Kfactors::K4fermion_normalized(mtsquared,TC);// K4_fermion_interp(mtsquared,TC,PyLoad)/norm1;
    double K5t = Kfactors::K5fermion_normalized(mtsquared,TC);// K5_fermion_interp(mtsquared,TC,PyLoad)/norm2;
    double K6t = Kfactors::K6fermion_normalized(mtsquared,TC);// K6_fermion_interp(mtsquared,TC,PyLoad)/norm2;
    double K8t = Kfactors::K8fermion_normalized(mtsquared,TC);// K8_fermion_interp(mtsquared,TC,PyLoad)/norm1;
    double K9t = Kfactors::K9fermion_normalized(mtsquared,TC);// K9_fermion_interp(mtsquared,TC,PyLoad)/norm1;

    double K1b = Kfactors::K1fermion_normalized(0,TC);// K1_fermion_interp(0,TC,PyLoad)/norm1;
    double K4b = Kfactors::K4fermion_normalized(0,TC);// K4_fermion_interp(0,TC,PyLoad)/norm1;
    double K5b = Kfactors::K5fermion_normalized(0,TC);// K5_fermion_interp(0,TC,PyLoad)/norm2;

    double K1h = Kfactors::K1boson_normalized(0,TC);// K1_boson_interp(0,TC,PyLoad)/norm1;
    double K4h = Kfactors::K4boson_normalized(0,TC);// K4_boson_interp(0,TC,PyLoad)/norm1;
    double K5h = Kfactors::K5boson_normalized(0,TC);// K5_boson_interp(0,TC,PyLoad)/norm2;


    double Dh = 20.0/TC;
    double Dt = 6.0/TC;

    double GH = mws/(50.0*TC);
    double GY = 4.2e-3*TC;
    double GM = mtsquared/(63*TC);
    double GSS = 4.9e-4*TC;



    double GHTot = K4h/(K1h*Dh);
    double GTTot = K4t/(K1t*Dt);
    double GBTot = K4b/(K1b*Dt);
    double GW = GHTot;


    double dtheta = calculate_theta(z,1);
    double d2theta = calculate_theta(z,2);

    double St = -vw * K8t * (mtsquared * d2theta + dmtsquared * dtheta) + vw * K9t * dtheta * mtsquared * dmtsquared;


    if(UseVelocityTransportEquations){
        double mut2 = x[0];
        double mub2 = x[1];
        double mutc2 = x[2];
        double muh2 = x[3];
        double ut2 = x[4];
        double ub2 = x[5];
        double utc2 = x[6];
        double uh2 = x[7];


        dxdt[0] = (double) ((-3 * K2t * K5t * mut2 * vw * vw * dmtsquared + 27 * GSS * K1b * K5t * mub2 * vw + 27 * GSS * K1t * K5t * mut2 * vw - 27 * GSS * K1t * K5t * mutc2 * vw + 6 * GM * K5t * mut2 * vw + 6 * GM * K5t * mutc2 * vw + 3 * GSS * K5t * mub2 * vw + 3 * GSS * K5t * mut2 * vw + 3 * GSS * K5t * mutc2 * vw - 3 * GW * K5t * mub2 * vw + 3 * GW * K5t * mut2 * vw + 3 * GY * K5t * muh2 * vw + 3 * GY * K5t * mut2 * vw + 3 * GY * K5t * mutc2 * vw + 3 * K6t * ut2 * vw * dmtsquared + 3 * GTTot * ut2 - St) / (K1t * K5t * vw * vw + K4t)) / 0.3e1;
        dxdt[1] = (9 * GSS * K1b * K5b * mub2 * vw + 9 * GSS * K1t * K5b * mut2 * vw - 9 * GSS * K1t * K5b * mutc2 * vw + GSS * K5b * mub2 * vw + GSS * K5b * mut2 * vw + GSS * K5b * mutc2 * vw + GW * K5b * mub2 * vw - GW * K5b * mut2 * vw + GY * K5b * mub2 * vw + GY * K5b * muh2 * vw + GY * K5b * mutc2 * vw + GBTot * ub2) / (K1b * K5b * vw * vw + K4b);
        dxdt[2] = (double) ((-3 * K2t * K5t * mutc2 * vw * vw * dmtsquared + 27 * GSS * K1b * K5t * mub2 * vw + 27 * GSS * K1t * K5t * mut2 * vw - 27 * GSS * K1t * K5t * mutc2 * vw + 6 * GM * K5t * mut2 * vw + 6 * GM * K5t * mutc2 * vw + 3 * GSS * K5t * mub2 * vw + 3 * GSS * K5t * mut2 * vw + 3 * GSS * K5t * mutc2 * vw + 3 * GY * K5t * mub2 * vw + 6 * GY * K5t * muh2 * vw + 3 * GY * K5t * mut2 * vw + 6 * GY * K5t * mutc2 * vw + 3 * K6t * utc2 * vw * dmtsquared + 3 * GTTot * utc2 - St) / (K1t * K5t * vw * vw + K4t)) / 0.3e1;
        dxdt[3] = (double) ((4 * GH * K5h * muh2 * vw + 3 * GY * K5h * mub2 * vw + 6 * GY * K5h * muh2 * vw + 3 * GY * K5h * mut2 * vw + 6 * GY * K5h * mutc2 * vw + 4 * GHTot * uh2) / (K1h * K5h * vw * vw + K4h)) / 0.4e1;
        dxdt[4] = (double) ((-3 * K1t * K6t * ut2 * vw * vw * dmtsquared - 3 * K2t * K4t * mut2 * vw * dmtsquared + 27 * GSS * K1b * K4t * mub2 + 27 * GSS * K1t * K4t * mut2 - 27 * GSS * K1t * K4t * mutc2 - 3 * GTTot * K1t * ut2 * vw + 6 * GM * K4t * mut2 + 6 * GM * K4t * mutc2 + 3 * GSS * K4t * mub2 + 3 * GSS * K4t * mut2 + 3 * GSS * K4t * mutc2 - 3 * GW * K4t * mub2 + 3 * GW * K4t * mut2 + 3 * GY * K4t * muh2 + 3 * GY * K4t * mut2 + 3 * GY * K4t * mutc2 + K1t * St * vw) / (K1t * K5t * vw * vw + K4t)) / 0.3e1;
        dxdt[5] = -(GBTot * K1b * ub2 * vw - 9 * GSS * K1b * K4b * mub2 - 9 * GSS * K1t * K4b * mut2 + 9 * GSS * K1t * K4b * mutc2 - GSS * K4b * mub2 - GSS * K4b * mut2 - GSS * K4b * mutc2 - GW * K4b * mub2 + GW * K4b * mut2 - GY * K4b * mub2 - GY * K4b * muh2 - GY * K4b * mutc2) / (K1b * K5b * vw * vw + K4b);
        dxdt[6] = (double) ((-3 * K1t * K6t * utc2 * vw * vw * dmtsquared - 3 * K2t * K4t * mutc2 * vw * dmtsquared + 27 * GSS * K1b * K4t * mub2 + 27 * GSS * K1t * K4t * mut2 - 27 * GSS * K1t * K4t * mutc2 - 3 * GTTot * K1t * utc2 * vw + 6 * GM * K4t * mut2 + 6 * GM * K4t * mutc2 + 3 * GSS * K4t * mub2 + 3 * GSS * K4t * mut2 + 3 * GSS * K4t * mutc2 + 3 * GY * K4t * mub2 + 6 * GY * K4t * muh2 + 3 * GY * K4t * mut2 + 6 * GY * K4t * mutc2 + K1t * St * vw) / (K1t * K5t * vw * vw + K4t)) / 0.3e1;
        dxdt[7] = (double) ((-4 * GHTot * K1h * uh2 * vw + 4 * GH * K4h * muh2 + 3 * GY * K4h * mub2 + 6 * GY * K4h * muh2 + 3 * GY * K4h * mut2 + 6 * GY * K4h * mutc2) / (K1h * K5h * vw * vw + K4h)) / 0.4e1;
    }
    else{
        double mut2 = x[0];
        double mub2 = x[1];
        double mutc2 = x[2];
        double muh2 = x[3];
        double dmut2 = x[4];
        double dmub2 = x[5];
        double dmutc2 = x[6];
        double dmuh2 = x[7];

        dxdt[0] = dmut2;
        dxdt[1] = dmub2;
        dxdt[2] = dmutc2;
        dxdt[3] = dmuh2;

        double dSt  = -vw*K8t *( 2*dmtsquared * d2theta) + vw*K9t * (d2theta *mtsquared + dtheta *dmtsquared)*dmtsquared;

        double ddmut2,ddmutc2,ddmuh2,ddmub2;

        ddmut2 =  ((-3 * K2t * K6t * mut2 * vw * vw * dmtsquared * dmtsquared + 27 * GSS * K1b * K6t * mub2 * vw * dmtsquared + 27 * GSS * K1t * K6t * mut2 * vw * dmtsquared - 27 * GSS * K1t * K6t * mutc2 * vw * dmtsquared - 3 * K1t * K6t * dmut2 * vw * vw * dmtsquared + 6 * GM * K6t * mut2 * vw * dmtsquared + 6 * GM * K6t * mutc2 * vw * dmtsquared + 3 * GSS * K6t * mub2 * vw * dmtsquared + 3 * GSS * K6t * mut2 * vw * dmtsquared + 3 * GSS * K6t * mutc2 * vw * dmtsquared - 3 * GTTot * K2t * mut2 * vw * dmtsquared - 3 * GW * K6t * mub2 * vw * dmtsquared + 3 * GW * K6t * mut2 * vw * dmtsquared + 3 * GY * K6t * muh2 * vw * dmtsquared + 3 * GY * K6t * mut2 * vw * dmtsquared + 3 * GY * K6t * mutc2 * vw * dmtsquared + 27 * GSS * GTTot * K1b * mub2 + 27 * GSS * GTTot * K1t * mut2 - 27 * GSS * GTTot * K1t * mutc2 - 3 * GTTot * K1t * dmut2 * vw + 6 * GM * GTTot * mut2 + 6 * GM * GTTot * mutc2 + 3 * GSS * GTTot * mub2 + 3 * GSS * GTTot * mut2 + 3 * GSS * GTTot * mutc2 - 3 * GTTot * GW * mub2 + 3 * GTTot * GW * mut2 + 3 * GTTot * GY * muh2 + 3 * GTTot * GY * mut2 + 3 * GTTot * GY * mutc2 - dSt) / K4t) / 0.3e1;
        ddmutc2 = ((-3 * K2t * K6t * mutc2 * vw * vw * dmtsquared * dmtsquared + 27 * GSS * K1b * K6t * mub2 * vw * dmtsquared + 27 * GSS * K1t * K6t * mut2 * vw * dmtsquared - 27 * GSS * K1t * K6t * mutc2 * vw * dmtsquared - 3 * K1t * K6t * dmutc2 * vw * vw * dmtsquared + 6 * GM * K6t * mut2 * vw * dmtsquared + 6 * GM * K6t * mutc2 * vw * dmtsquared + 3 * GSS * K6t * mub2 * vw * dmtsquared + 3 * GSS * K6t * mut2 * vw * dmtsquared + 3 * GSS * K6t * mutc2 * vw * dmtsquared - 3 * GTTot * K2t * mutc2 * vw * dmtsquared + 3 * GY * K6t * mub2 * vw * dmtsquared + 6 * GY * K6t * muh2 * vw * dmtsquared + 3 * GY * K6t * mut2 * vw * dmtsquared + 6 * GY * K6t * mutc2 * vw * dmtsquared + 27 * GSS * GTTot * K1b * mub2 + 27 * GSS * GTTot * K1t * mut2 - 27 * GSS * GTTot * K1t * mutc2 - 3 * GTTot * K1t * dmutc2 * vw + 6 * GM * GTTot * mut2 + 6 * GM * GTTot * mutc2 + 3 * GSS * GTTot * mub2 + 3 * GSS * GTTot * mut2 + 3 * GSS * GTTot * mutc2 + 3 * GTTot * GY * mub2 + 6 * GTTot * GY * muh2 + 3 * GTTot * GY * mut2 + 6 * GTTot * GY * mutc2 - dSt) / K4t) / 0.3e1;

        ddmuh2 = (GHTot * (-4 * vw * K1h * dmuh2 + 4 * GH * muh2 + 3 * mub2 * GY + 6 * GY * muh2 + 3 * mut2 * GY + 6 * GY * mutc2) / K4h) / 0.4e1;


        ddmub2 = GBTot * (9 * GSS * K1b * mub2 + 9 * GSS * K1t * mut2 - 9 * GSS * K1t * mutc2 - vw * K1b * dmub2 + GSS * mub2 + GSS * mut2 + GSS * mutc2 + mub2 * GW - mut2 * GW + GY * mub2 + muh2 * GY + mutc2 * GY) / K4b;


        dxdt[4] = ddmut2;
        dxdt[6] = ddmutc2;
        dxdt[7] = ddmuh2;
        dxdt[5] = ddmub2;

        //		  //ddmut2
        //		dxdt[4] = (double) ((-3 * K2t * K6t * mut2 * vw * vw * dmtsquared * dmtsquared + 27 * GSS * K1b * K6t * mub2 * vw * dmtsquared + 27 * GSS * K1t * K6t * mut2 * vw * dmtsquared - 27 * GSS * K1t * K6t * mutc2 * vw * dmtsquared - 3 * K1t * K6t * dmut2 * vw * vw * dmtsquared + 6 * GM * K6t * mut2 * vw * dmtsquared + 6 * GM * K6t * mutc2 * vw * dmtsquared + 3 * GSS * K6t * mub2 * vw * dmtsquared + 3 * GSS * K6t * mut2 * vw * dmtsquared + 3 * GSS * K6t * mutc2 * vw * dmtsquared - 3 * GTTot * K2t * mut2 * vw * dmtsquared - 3 * GW * K6t * mub2 * vw * dmtsquared + 3 * GW * K6t * mut2 * vw * dmtsquared + 3 * GY * K6t * muh2 * vw * dmtsquared + 3 * GY * K6t * mut2 * vw * dmtsquared + 3 * GY * K6t * mutc2 * vw * dmtsquared + 27 * GSS * GTTot * K1b * mub2 + 27 * GSS * GTTot * K1t * mut2 - 27 * GSS * GTTot * K1t * mutc2 - 3 * GTTot * K1t * dmut2 * vw + 6 * GM * GTTot * mut2 + 6 * GM * GTTot * mutc2 + 3 * GSS * GTTot * mub2 + 3 * GSS * GTTot * mut2 + 3 * GSS * GTTot * mutc2 - 3 * GTTot * GW * mub2 + 3 * GTTot * GW * mut2 + 3 * GTTot * GY * muh2 + 3 * GTTot * GY * mut2 + 3 * GTTot * GY * mutc2 - dSt) / K4t) / 0.3e1;
        //		//ddmutc2
        //		dxdt[6] = (double) ((-3 * vw * vw * K2t * dmtsquared * dmtsquared * mutc2 * K6t + 27 * GSS * K1b * K6t * mub2 * vw * dmtsquared + 27 * GSS * K1t * K6t * mut2 * vw * dmtsquared - 27 * GSS * K1t * K6t * mutc2 * vw * dmtsquared - 3 * vw * vw * K1t * dmutc2 * K6t * dmtsquared + 6 * GM * K6t * mut2 * vw * dmtsquared + 6 * GM * K6t * mutc2 * vw * dmtsquared + 3 * GSS * K6t * mub2 * vw * dmtsquared + 3 * GSS * K6t * mut2 * vw * dmtsquared + 3 * GSS * K6t * mutc2 * vw * dmtsquared + 3 * GY * K6t * vw * dmtsquared * mub2 + 6 * GY * K6t * muh2 * vw * dmtsquared + 3 * GY * K6t * mut2 * vw * dmtsquared + 6 * GY * K6t * mutc2 * vw * dmtsquared - dSt) / K4t) / 0.3e1;
        //		//ddmuh2
        //		dxdt[7] = (double) (GHTot * (-4 * vw * K1h * dmuh2 + 3 * mub2 * GY + 6 * muh2 * GY + 3 * mut2 * GY + 6 * mutc2 * GY + 4 * GH * muh2) / K4h) / 0.4e1;
        //		//ddmub2
        //		dxdt[5] = GBTot * (9 * GSS * mub2 * K1b + 9 * mut2 * GSS * K1t - 9 * GSS * mutc2 * K1t - vw * K1b * dmub2 + GSS * mub2 + mut2 * GSS + GSS * mutc2 + mub2 * GW - mut2 * GW + mub2 * GY + muh2 * GY + mutc2 * GY) / K4b;


    }

    return *this;

}


std::vector<double> calculateTransportEquation(
        const double& z,
        const std::vector<double>& parStart,
        const struct GSL_integration_mubl& params)
{
    using namespace boost::numeric::odeint;
    const double C_AbsErr = 1e-10; //1.0e-10
    const double C_RelErr = 1e-3; // 1.0e-6

    std::size_t dim = 8;
    // dim = 1;
    state_type x(dim);
    for(std::size_t i=0;i<dim;i++) x[i] = parStart[i];

    transport_equations transport(params);
    std::vector<state_type> x_vec;
    std::vector<double> times;

    double zInitial=params.getZMAX();

    double stepsize_initial=(z-zInitial)*1e-5;
    if(stepsize_initial == 0) stepsize_initial=-1e-7;
    double abs_err = C_AbsErr;
    double rel_err = C_RelErr;



    // steps = integrate(transport,x,zInitial,z,stepsize_initial);
    //	steps = integrate_adaptive( make_controlled( abs_err , rel_err , error_stepper_type() ) ,RGE , x , EnergyStart , EnergyEnd , stepsize_initial );

    integrate_adaptive( make_controlled( abs_err , rel_err , error_stepper_type() ) ,transport ,
                        x , zInitial , z , stepsize_initial,push_back_state_and_time( x_vec , times ) );

    std::vector<double> parEnd;
    for(std::size_t i=0;i<dim;i++) parEnd.push_back(x[i]);
    for(std::size_t i=0;i<dim;i++){
        if(std::abs(parEnd[i]) <= std::pow(10,-16)) parEnd[i] = 0;
    }

    std::vector<double> topres,vev;
    vev=transport.calculate_vev(z);
    topres=transport.get_top_mass_and_derivative(vev);
    double mtsquared = topres.at(0);
    double K1t = Kfactors::K1fermion_normalized(mtsquared,params.getTC());
    double K1b = Kfactors::K1fermion_normalized(0,params.getTC());

    parEnd.push_back(K1t);
    parEnd.push_back(K1b);

    if(std::isnan(parEnd.at(0))){
        std::cout << "Nan in " << __func__ << std::endl
                  << "parEnd.size() = " << parEnd.size() << "\nparEnd = ";
        for(std::size_t i=0;i<parEnd.size();i++) std::cout << parEnd.at(i) << sep;
        std::cout << std::endl;
    }
    return parEnd;
}



double mubl_func(double z,  void *p){
    struct GSL_integration_mubl * params = (struct GSL_integration_mubl *) p;
    std::vector<double> parStart,parEnd;
    for(int i=0;i<8;i++) parStart.push_back(0);

    parEnd=calculateTransportEquation(z,parStart,*params);
    double K1t = parEnd.at(8);
    double K1b = parEnd.at(9);
    double mut2 = parEnd.at(0);
    double mub2 = parEnd.at(1);
    double mutc2 = parEnd.at(2);

    double res = 0.5*(1+4*K1t)*mut2 + 0.5*(1+4*K1b)*mub2-2*K1t*mutc2;
    if(std::isnan(res)){
        std::cout << "res = nan at z = " << z << std::endl;
        std::cout << "K1t = " << K1t << "\nK1b = " << K1b << "\nmut2 = " << mut2 << "\nmub2 = " << mub2
                  << "\nmutc2 = " << mutc2 << std::endl;
    }
    return res;
}



double eta_integrand_func(double z,  void *p){
    struct GSL_integration_mubl * params = static_cast<GSL_integration_mubl*>(p);
    double mubl = mubl_func(z,p);
    double GWS = 1e-6*params->getTC();
    double nu = 45*GWS/(4*params->getvw());
    double damping = std::exp(-nu*z);
    double res = mubl*damping;
    return res;
}


double Integrate_mubl(const struct GSL_integration_mubl& pIn){
    std::size_t workspace_size = 1000;
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(workspace_size);
    double result, error;
    GSL_integration_mubl p = pIn;
    gsl_function F;
    F.function=&eta_integrand_func;
    F.params=&p;

    double zmax = p.getZMAX();

    gsl_integration_qags(&F,0,zmax,0,1e-7,workspace_size,w,&result,&error);

    gsl_integration_workspace_free(w);


    return result;
}


GSL_mubl_interpolation generate_mubl_spline(const struct GSL_integration_mubl& p, int nstep)
{
    double MaxZ = p.getZMAX();
    double stepsize= (MaxZ-0)/nstep;
    std::vector<double> ydata;
    struct GSL_integration_mubl pIn;

    pIn = p;
    for(int i=0;i<=nstep;i++){
        double z = i*stepsize;
        double resmubl = mubl_func(z,&pIn);
        ydata.push_back(resmubl);

    }

    boost::math::cubic_b_spline<double> splinef(ydata.data(),ydata.size(),0,stepsize);
    GSL_mubl_interpolation spline;
    spline.spline = splinef;
    spline.vw = p.getvw();
    spline.TC = p.getTC();
    return spline;
}


double mubl_interpolation(double z,void *p){
    struct GSL_mubl_interpolation * params = static_cast<GSL_mubl_interpolation*>(p);
    double GWS = 1e-6*params->TC;
    double nu = 45*GWS/(4*params->vw);
    double res = params->spline(z)*std::exp(-nu*z);
    return res;
}


double Integrate_mubl_interpolated(const struct GSL_integration_mubl& p){

    double vw = p.getvw();

    int nstep=100;
    auto spline = generate_mubl_spline(p, nstep);

    std::size_t workspace_size = 1000;
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(workspace_size);
    double result, error;
    gsl_function F;
    F.function=&mubl_interpolation;
    F.params=&spline;

    double zmax = p.getZMAX();

    gsl_integration_qags(&F,0,zmax,0,1e-9,workspace_size,w,&result,&error);

    gsl_integration_workspace_free(w);


    double GWS = 1e-6*p.getTC();
    double gstar = 106.75;
    double prefac = 405*GWS/(4*M_PI*M_PI*vw*gstar*p.getTC());

    result *= prefac;
    return result;
}

}
}
