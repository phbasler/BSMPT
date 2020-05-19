/*
 * gen_func_fluid.cpp
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

#include <BSMPT/baryo_calculation/Fluid_Type/gen_func_fluid.h>
#include <BSMPT/models/ClassPotentialOrigin.h>
#include <BSMPT/utility.h>

/**
 * @file
 */

namespace BSMPT{
namespace Baryo{
/**
 * @brief C_AbsErr Absolute tolerance for boost::integrate_adaptive
 */
const double C_AbsErr = 1e-6;
/**
 * @brief C_RelErr Relative tolerance for boost::integrate_adaptive
 */
const double C_RelErr = 1e-9;

void Calc_Gam_M::operator()(const state_type &Gam, state_type &dGam,const double k)
{
    (void) Gam;
    bool debug = false;
    double fsqrt = 2 *Calc_Gam_M::m_full*Calc_Gam_M::m_full;
    double Nc = 3;
    double pi = M_PI;
    double w = std::sqrt(std::pow(k,2)+ msqrt_thermal);
    double delta_w = dmsqrt_thermal/(2*std::sqrt(std::pow(k,2)+msqrt_thermal));
//    double delta_w = 0;


    double integrand = 0 ;
    if(std::abs(msqrt_thermal)>C_threshold)
    {
        integrand+=std::pow(k,2)/GamT;
        integrand-=std::pow(k,4)/(std::pow(w,2)*GamT);
        integrand+=std::pow(k,4)*GamT/(std::pow(w,4));
        integrand+=delta_w*std::pow(k,4)/(std::pow(w,3)*GamT);
        integrand-=2*delta_w*std::pow(k,4)*GamT/std::pow(w,5);
        integrand=integrand*hf(w);
        //Checked with Mathematica Calc_Relax.m-->GammaMinus
    }
    if(std::abs(msqrt_thermal) < C_threshold)
    {
        throw std::runtime_error("Thermal mass under numerical threshold-->Proper limit needed in Calc_Gam_M");
        //Check with Mathematica Calc_Relax.m-->GammaMinus
    }
    if(std::abs(integrand)>C_smallcut){
        dGam[0] = (6./std::pow(T,2))*(Nc*fsqrt * integrand)/(2*pi*pi*T) ;
    }
    if(std::abs(integrand)<C_smallcut) dGam[0]=0;
    if(debug)
    {
        std::cout<<"fsqrt = "<<fsqrt<<std::endl;
        std::cout<<"T = "<<T<<std::endl;
        std::cout<<"beta = "<<beta<<std::endl;
        std::cout<<"GamT = "<<GamT<<std::endl;
        std::cout<<"w = "<<w<<std::endl;
        std::cout<<"k = "<<k<<std::endl;
        std::cout<<"dGam[0] = "<<dGam[0]<<std::endl;
        std::cout<<"integrand = " <<integrand<<std::endl;
        std::cout<<"\n\n";
    }
}

void Calc_Gam_M::set_class(double T_in, double vw_in, double m_in, double msqrt_thermal_in, double dmsqrt_thermal_in)
{
    bool debug =  false;
    if(debug) std::cout<<"Begin of debug in "<<__func__<<std::endl;
    Calc_Gam_M::T = T_in;
    Calc_Gam_M::vw= vw_in;
    Calc_Gam_M::m_full=m_in;
    Calc_Gam_M::beta=1/T;
    Calc_Gam_M::GamT=0.16*T;
    Calc_Gam_M::msqrt_thermal = msqrt_thermal_in;
    Calc_Gam_M::dmsqrt_thermal = dmsqrt_thermal_in;

    if(debug)
    {
        std::cout<<"set parameters in Calc_Gam_M:\n";
        std::cout<< sep << "T -->"<<Calc_Gam_M::T<<std::endl;
        std::cout<< sep << "vw -->"<<Calc_Gam_M::vw<<std::endl;
        std::cout<< sep << "m_full -->"<<Calc_Gam_M::m_full<<std::endl;
        std::cout<< sep << "m_thermal-->"<<Calc_Gam_M::msqrt_thermal<<std::endl;
        std::cout<< sep << "dm_thermal-->"<<Calc_Gam_M::dmsqrt_thermal<<std::endl;
        std::cout<< sep << "beta --> "<<Calc_Gam_M::beta<<std::endl;
        std::cout<< sep << "GamT --> "<<Calc_Gam_M::GamT<<std::endl;
    }
}

double Calc_Gam_M::hf(double w)
{
    double beta_calc = Calc_Gam_M::beta;
    double exp = std::exp(w*beta_calc);
    if(std::isnan(exp)) return 0;
    double res = exp / std::pow(1+exp,2);
    if( std::abs(res) < C_smallcut) res = 0;
    if(std::isnan(res)) throw std::runtime_error("NaN @ Calc_Gam_M::hf");
    return res;
}
double Calc_Gam_M::hf_prime(double w)
{
    bool debug = false;
    double beta_calc = Calc_Gam_M::beta;   
    double exp = std::exp(w*beta_calc);
    double res = - beta*exp*(exp-1)/std::pow((1+exp),3);
    if(debug) std::cout<<sep << sep <<" exp0 = " << res<<std::endl;
    if(std::abs(res) < C_smallcut) res =0;
    if(debug) std::cout<<sep << sep <<" exp1 = " << res<<std::endl;
    return res;
}

double Nintegrate_GamM(Calc_Gam_M& C_Gam)
{
    bool debug = false ;
    if(debug) std::cout<<"Begin of debug in "<<__func__<<std::endl;
    state_type Gam(1);
    //Paramters for the numerical evaluation
        double stepsize_initial = 0.01;
        double k_start = 0;
        //upper limit on k is chosen such that no std::limit overflow occurs in hf, hf_prime functions
        double k_end = 1e4;
    if(debug){
        std::cout<<"\tCalled "<<__func__<<" with the following arguments:\n";
        std::cout<<"\t\tm = "<<C_Gam.m_full<<std::endl;
        std::cout<<"\t\tm_thermal = "<<C_Gam.msqrt_thermal<<std::endl;
        std::cout<<"\t\tdm_thermal = "<<C_Gam.dmsqrt_thermal<<std::endl;
        std::cout<<"\t\tT = "<<C_Gam.T<<std::endl;
        std::cout<<"\t\tvw = "<<C_Gam.vw<<std::endl;
    }
    double abs_err = C_AbsErr;
    double rel_err = C_RelErr;
    // typedef runge_kutta_cash_karp54< state_type > error_stepper_type;
    integrate_adaptive(make_controlled( abs_err , rel_err , error_stepper_type() ) , C_Gam , Gam , k_start , k_end , stepsize_initial );
    if(debug) std::cout<<"\tGam_M = " << Gam[0]<<std::endl;
    if(debug) std::cout<<"End: "<<__func__<<std::endl;
    if(std::isnan(Gam[0]))
    {
        throw std::runtime_error("NaN in Nintegrate_GamM"); 
    }
    else return Gam[0];
}


void Calc_Scp::operator()(const state_type &Scp , state_type &dScp , const double k)
{
    bool debug = false;
    state_type res = Scp;
    double pi = M_PI;
    double NC = 3; 
    double w = std::sqrt(std::pow(k,2)+msqrt_thermal);
    double delta_w = dmsqrt_thermal/(2*std::sqrt(k*k+msqrt_thermal));
    double GamT= 0.16*T;
    double Imff = 2 * Calc_Scp::m_full*Calc_Scp::m_full*theta_prime;
    
    double integrand = 0;
    if(std::abs(Calc_Scp::msqrt_thermal) > 1e-3)
    {
        integrand-=std::pow(k,4)*GamT/(2*std::pow(w,5));
        if(std::abs(integrand)<C_smallcut)integrand=0;
        integrand+=5*std::pow(k,4)*GamT*delta_w/(4*std::pow(w,6));
        if(std::abs(integrand)<C_smallcut)integrand=0;
        integrand+=std::pow(k,4)*GamT*nf(w)/std::pow(w,5);
        if(std::abs(integrand)<C_smallcut)integrand=0;
        integrand-=5*std::pow(k,4)*GamT*delta_w*nf(w)/(2*std::pow(w,6));
        if(std::abs(integrand)<C_smallcut)integrand=0;
        integrand+= (k*k/(2*GamT) - std::pow(k,4)/(2*w*w*GamT) - std::pow(k,4)*GamT/(2*std::pow(k,4)))*nf_prime(w);
        if(std::abs(integrand)<C_smallcut)integrand=0;
        integrand+=(std::pow(k,4)/(2*std::pow(w,3)*GamT)+(3*std::pow(k,4)*GamT)/(2*std::pow(w,5)))*delta_w*nf_prime(w);
        if(std::abs(integrand)<C_smallcut)integrand=0;
        //Checked with Mathematica Calc_Relax.m-->Scp        
    }
    if(std::abs(Calc_Scp::msqrt_thermal) < 1e-3)
    {
        throw std::runtime_error("No valid thermal mass-->proper limit needed in Calc_Scp");
    }
    if(std::abs(integrand)>C_smallcut){
         dScp[0] = (NC*vw*Imff*integrand)/(pi*pi) ;
    }
    else dScp[0]=0;

    if(debug)
    {
        std::cout<<sep << "m = "<<m_full<<std::endl;
        std::cout<<sep << "m_thermal = " <<msqrt_thermal<<std::endl;
        std::cout<<sep << "dm_thermal = " <<dmsqrt_thermal<<std::endl;
        std::cout<<sep << "T = "<<T<<std::endl;
        std::cout<<sep << "vw = "<<vw<<std::endl;
        std::cout<<sep << "theta_prime  = "<<theta_prime<<std::endl;
        std::cout<<sep << "Imff = "<<Imff<<std::endl;
        std::cout<<sep << "integrand = "<<integrand<<std::endl;
        std::cout<<sep << "dScp[0] = "<<dScp[0]<<std::endl;
    }
}

double Calc_Scp::nf(double w)
{
    double res = 1./(std::exp(w/Calc_Scp::T)+1);
    if(std::abs(res)<C_smallcut) res = 0;
    return res;
}
double Calc_Scp::nf_prime(double w)
{
    double beta = 1./Calc_Scp::T;
    double res =  -(beta * std::exp(beta*w))/std::pow(1+std::exp(beta*w),2);
    if(std::abs(res)<C_smallcut) res = 0; 
    return res;
}
void Calc_Scp::set_class(double T_in, double vw_in, double mt_in ,double theta_prime_in , double msqrt_thermal_in,double dmsqrt_thermal_in)
{
    Calc_Scp::m_full = mt_in;
    Calc_Scp::vw = vw_in;
    Calc_Scp::theta_prime = theta_prime_in;
    Calc_Scp::T = T_in;
    Calc_Scp::msqrt_thermal=msqrt_thermal_in;
    Calc_Scp::dmsqrt_thermal=dmsqrt_thermal_in;
}

double Nintegrate_Scp(Calc_Scp& C_Scp)
{
    bool debug = false;
    if(debug) std::cout<<"Begin of debug in "<<__func__<<std::endl;
    // typedef runge_kutta_cash_karp54< state_type > error_stepper_type;
	// typedef controlled_runge_kutta< error_stepper_type > controlled_stepper_type;
    state_type Scp(1);

    double stepsize_initial = 0.01;
    double abs_err = C_AbsErr;
    double rel_err = C_RelErr;
    double k_start = 0;
    if(debug){
        std::cout<<"\tCalled "<<__func__<<"with the following arguments:\n";
        std::cout<<"\t\tm = "<<C_Scp.m_full<<std::endl;
        std::cout<<"\t\tT = "<<C_Scp.T<<std::endl;
        std::cout<<"\t\tvw = "<<C_Scp.T<<std::endl;
        std::cout<<"\t\ttheta_prime = " <<C_Scp.theta_prime<<std::endl;
        std::cout<<"END: "<<__func__<<std::endl;
    }

    //Ranges checked with NIntegrage from mathematica 
    double k_end = 1e4;
    if(debug) std::cout<<"Begin of numerical integration"<<std::endl;
    integrate_adaptive(make_controlled( abs_err , rel_err , error_stepper_type() ) , C_Scp , Scp , k_start , k_end , stepsize_initial );
	if(debug) std::cout<<"numerical evaluation Scp = " << Scp[0]<<std::endl;

    if(std::isnan(Scp[0]))
    {
        throw std::runtime_error("NaN in Nintegrate_SCP");
    }
    else return Scp[0];

}

void Calc_kappa_t::set_class(double Temp_in,double mt_in)
{
    Calc_kappa_t::Temp = Temp_in;
    Calc_kappa_t::mt = mt_in;
}
void Calc_kappa_t::operator()(const state_type& kappa, state_type& dkappa,const double x)
{
    //To avoid compiler errors
    state_type res  = kappa;
    double cF = 6;
    double prefactor = cF/(M_PI*M_PI);
    dkappa[0] = prefactor*(x*std::exp(x)*std::sqrt(x*x-std::pow(mt/Temp,2))/std::pow(std::exp(x)+1,2)) ;
}

double NIntegrate_kappa(const Calc_kappa_t& C_kap)
{
    bool debug = false;
    if(debug) std::cout<<"Begin of debug in "<<__func__<<std::endl;

    state_type kappa(1);

    double stepsize_initial = 1e-5;
    double abs_err = C_AbsErr;
    double rel_err = C_RelErr;
    double eps = 1e-5;
    double xstart = C_kap.mt/C_kap.Temp+eps;
    double xend = 3*1e2;

    integrate_adaptive(make_controlled( abs_err , rel_err , error_stepper_type() ) , C_kap , kappa , xstart , xend , stepsize_initial );
	if(debug) std::cout<<"numerical evaluation kappa = " << kappa[0]<<std::endl;

    if(std::isnan(kappa[0])) throw std::runtime_error("NaN in NIntegrate_kappa");
    return kappa[0];
  
}

void Calc_eta::operator()(const state_type& eta, state_type& deta, const double z)
{
    //avoid compiler errors 
    (void) eta;
    //prefactor is not used here for numerical stability
    deta[0] =  Calc_eta::nL_cub(z) * std::exp(exponent_prefactor*z);

}
void Calc_eta::set_class(std::vector<double> array_z_in,std::vector<double> array_nL_in , double Temp_in,double vw_in)
{
    bool debug = false;
	if(debug) std::cout<<"Begin of debug in "<<__func__<<std::endl;

    Calc_eta::array_z = array_z_in;
    Calc_eta::array_nL = array_nL_in;
    Calc_eta::Temp = Temp_in;
    Calc_eta::vw = vw_in;

    double R =15./4;
    double kappa = 20;
    double gBMP=0.65;
    double gstar=106.75;
    double s_entr=gstar*std::pow(Temp,3)*2*M_PI*M_PI/(45) ;
    double Dq = 6./Temp;
    double alpha_W = gBMP*gBMP/(4*M_PI);
    double GamWS = 6*kappa*std::pow(alpha_W,5)*Temp;
    double alpha_P = (vw+std::sqrt(4*Dq*GamWS*R+vw*vw))/(2*Dq);
    double alpha_M = (vw-std::sqrt(4*Dq*GamWS*R+vw*vw))/(2*Dq);
    double omega = std::sqrt(vw*vw + 4*GamWS*Dq*R);
    //FULL EXPRESSION (Calculated from scratch)
    prefactor = -3*GamWS/(2*omega*s_entr);
    exponent_prefactor = (vw-omega)/(2*Dq);
    // exponent_prefactor = -alpha_M;
    //SIMPLIFIED EXPRESSION
    // prefactor = -3*GamWS/(2*vw*s_entr);
    // exponent_prefactor = -R*GamWS/vw;
    double t0,h;
    if(array_z.at(0)<array_z.at(array_z.size()-1))
    {
        t0 = array_z.at(0);
        h = array_z.at(1)-array_z.at(0);
    }
    if(array_z.at(0)>array_z.at(array_z.size()-1))
    {
        t0 = array_z.at(array_z.size()-1);
        h = array_z.at(0)-array_z.at(1);
    }
	boost::math::cubic_b_spline<double> spline(array_nL.begin(), array_nL.end(), t0, h);
	Calc_eta::nL_cub = spline;
    if(debug)
    {
        std::cout<<"Ranges of the spline:\n";
        std::cout<<"z in ["<<array_z[0]<<","<<sep<<array_z[99]<<"]\n";
        std::cout<<"nL in ["<<array_nL[0]<<","<<sep<<array_nL[99]<<"]\n";
    }
    if(debug)
    {
        std::cout<<sep << "Temp = "<< Temp<<std::endl;
        std::cout<<sep << "vw = "<< vw<<std::endl;
        std::cout<<sep << "prefactor = "<< prefactor<<std::endl;
        std::cout<<sep << "exponent_prefactor = "<< exponent_prefactor<<std::endl;
        std::cout<<sep << "GamWS = "<<GamWS<<std::endl;
        std::cout<<sep << "s_entr = "<<s_entr<<std::endl;
        std::cout<<sep << "alpha_M = "<<alpha_M<<std::endl;
        std::cout<<sep << "alpha_P = "<<alpha_P<<std::endl;
        std::cout<<"test call spline nL(0) =" << nL_cub(array_z[0])<< sep << "exact: nL = "<<array_nL[0]<<std::endl;
        std::cout<<"test call spline nL(2) =" << nL_cub(array_z[2])<< sep << "exact: nL = "<<array_nL[2]<<std::endl;
        std::cout<<"test call spline nL(10) =" << nL_cub(array_z[10])<< sep << "exact: nL = "<<array_nL[10]<<std::endl;
    }
}
void Calc_eta::set_class(std::pair< std::vector<double> , std::vector<double>> arr_in, double Temp_in,double vw_in)
{
    bool debug = false;
	if(debug) std::cout<<"Begin of debug in "<<__func__<<std::endl;

    Calc_eta::array_z = arr_in.first;
    Calc_eta::array_nL = arr_in.second;
    Calc_eta::Temp = Temp_in;
    Calc_eta::vw = vw_in;

    double R =15./4;
    double kappa = 20;
    double gBMP=0.65;
    double gstar=106.75;
    double s_entr=gstar*std::pow(Temp,3)*2*M_PI*M_PI/(45.) ;
    double Dq = 6./Temp;
    double alpha_W = gBMP*gBMP/(4*M_PI);
    double GamWS = 6*kappa*std::pow(alpha_W,5)*Temp; 
    double alpha_P = (vw+std::sqrt(4*Dq*GamWS*R+vw*vw))/(2*Dq);
    double alpha_M = (vw-std::sqrt(4*Dq*GamWS*R+vw*vw))/(2*Dq);
    double omega = std::sqrt(vw*vw + 4*GamWS*Dq*R);
    //FULL EXPRESSION (Calculated from scratch)
    prefactor = -3*GamWS/(2*omega*s_entr);
    exponent_prefactor = (vw-omega)/(2*Dq);
    // exponent_prefactor = -alpha_M;
    //SIMPLIFIED EXPRESSION
    // prefactor = -3*GamWS/(2*vw*s_entr);
    // exponent_prefactor = -R*GamWS/vw;
    double t0,h;
    if(array_z.at(0)<array_z.at(array_z.size()-1))
    {
        t0 = array_z.at(0);
        h = array_z.at(1)-array_z.at(0);
    }
    if(array_z.at(0)>array_z.at(array_z.size()-1))
    {
        t0 = array_z.at(array_z.size()-1);
        h = array_z.at(0)-array_z.at(1);
    }
	boost::math::cubic_b_spline<double> spline(array_nL.begin(), array_nL.end(), t0, h);
	Calc_eta::nL_cub = spline;

    if(debug)
    {
        std::cout<<sep << "Temp = "<< Temp<<std::endl;
        std::cout<<sep << "vw = "<< vw<<std::endl;
        std::cout<<sep << "prefactor = "<< prefactor<<std::endl;
        std::cout<<sep << "exponent_prefactor = "<< exponent_prefactor<<std::endl;
        std::cout<<sep << "GamWS = "<<GamWS<<std::endl;
        std::cout<<sep << "s_entr = "<<s_entr<<std::endl;
        std::cout<<sep << "alpha_M = "<<alpha_M<<std::endl;
        std::cout<<sep << "alpha_P = "<<alpha_P<<std::endl;
        std::cout<<"test call spline nL(0) =" << nL_cub(array_z[0])<< sep << "exact: nL = "<<array_nL[0]<<std::endl;
        std::cout<<"test call spline nL(2) =" << nL_cub(array_z[2])<< sep << "exact: nL = "<<array_nL[2]<<std::endl;
        std::cout<<"test call spline nL(10) =" << nL_cub(array_z[10])<< sep << "exact: nL = "<<array_nL[10]<<std::endl;
    }
}
double Nintegrate_eta(const Calc_eta& C_eta, const double &z_start, const double &z_end)
{
    bool debug = false;
    if(debug) std::cout<<"Begin of debug in "<<__func__<<std::endl;
    state_type eta(1);

    double stepsize_initial = -1e-5;
    double abs_err = C_AbsErr;
    double rel_err = C_RelErr;

    if(debug) std::cout<<"Begin of numerical integration"<<std::endl;
    integrate_adaptive( make_controlled( abs_err , rel_err , error_stepper_type() ) , C_eta , eta , z_end , z_start , stepsize_initial );
	if(debug) std::cout<<"numerical evaluation Scp = " << eta[0]<<std::endl;

    if( std::isnan(eta[0]) or std::isnan(C_eta.prefactor)) throw std::runtime_error("NaN in Nintegrate_eta");
    return C_eta.prefactor*eta[0];

}


void gen_fluid::top_func(double z , std::vector<double>& m_quark, std::vector<double>& m_quark_prime)
{
  bool debug = false;
  if(debug) std::cout<<"Begin debug in function " <<__func__<<std::endl;
  std::vector<double> m_i;
  std::vector<double> dif_mt;
  std::vector<double> dif_mb;
  dif_mt.clear();
  dif_mb.clear();
  std::size_t nquark = gen_fluid::modelPointer->get_NQuarks();
  std::vector<double> phi,gen_phi;
  phi.clear();
  gen_phi.clear();
  std::size_t gen_nvev = gen_fluid::gen_vcritical.size();
  for(size_t i =0 ; i<gen_nvev;i++) gen_phi.push_back(gen_fluid::gen_vcritical.at(i)/2*(1-std::tanh(z/gen_fluid::LW)));


  for(size_t i=1;i<=gen_nvev;i++)
    {
	  m_i.clear();
      m_i=gen_fluid::modelPointer->QuarkMassesSquared(gen_phi,i);
	  /*
	  m_i: nquark eigenvalues and nquark derivatives with respect to field i
	  top is the heaviest quark -> at position nquark-1 and the derivative at position 2 nquark -1
      bot is the second heaviest quark ->at position nquark-2 and the derivative at position 2 nquark -2 
	  */
	  dif_mt.push_back(m_i.at(2*nquark-1));
      //top --> anti_top-->b-->anti b 
      dif_mb.push_back(m_i.at(2*nquark-3));
    }
    double mtop = std::sqrt(std::abs(m_i.at(nquark-1)));
    double mbot = std::sqrt(std::abs(m_i.at(nquark-3)));
    double res_t=0;
    double res_b=0;
    std::vector<double> grad_omega = gen_fluid::omegaprime(z);
        if(debug) std::cout<<std::endl;
        if(debug) for(size_t i=0;i<gen_nvev;i++) std::cout<<"gen_vev_"<<i<<" = "<<gen_phi.at(i)<<std::endl;
        if(debug) for(size_t i=0;i<grad_omega.size();i++)std::cout<<"grad_omega_"<<i<<" = "<< grad_omega.at(i)<<"\n";
        if(debug) std::cout<<std::endl;
        if(debug) for(size_t i=0;i<dif_mt.size();i++)std::cout<<"dif_mt_"<<i<<" ="<< dif_mt.at(i)<<"\n";
        if(debug) std::cout<<std::endl;
    for(size_t i=0;i<dif_mt.size();i++) res_t+=dif_mt.at(i)*grad_omega.at(i);
    for(size_t i=0;i<dif_mb.size();i++) res_b+=dif_mb.at(i)*grad_omega.at(i);

    if(debug) std::cout<<"msqrtprime = "<<res_t<<std::endl;
    if(std::isnan(res_t)) throw std::runtime_error("NAN @ msqrtprime");
    if(std::isnan(res_b)) throw std::runtime_error("NAN @ msqrtprime");

    m_quark.resize(2);
    m_quark_prime.resize(2);
    
    m_quark[0] = mtop;
    m_quark[1] = mbot;
    m_quark_prime[0] = res_t;
    m_quark_prime[1] = res_b;
}
std::vector<double> gen_fluid::omegaprime(double z)
{
  std::vector<double> dummy;
  dummy.clear();
  for(size_t i=0;i<gen_fluid::gen_vcritical.size();i++)
  {
	dummy.push_back((-gen_fluid::gen_vcritical.at(i)/(2*gen_fluid::LW))*(1./(std::pow(std::cosh(z/gen_fluid::LW),2))));
  }
  return dummy;
}


double gen_fluid::atan2_mod(double Im , double Re)
{
  double foo = std::atan2(Im,Re);
  if(std::abs(foo)<M_PI/2.) return  std::atan2(Im,Re);
  else return  std::atan2(-Im, - Re);
}


//Calc_Theta with the generic CP-violating phase input 
std::vector<double> gen_fluid::Calc_theta(double z , double CP_sym , double CP_brk)
{
    std::vector<double> res;
    res.resize(3);
    if(modelPointer->get_Model() != ModelID::ModelIDs::C2HDM ) std::cerr << "This is only programmed for the C2HDM" << std::endl;
    bool debug = false;
    if(debug) std::cout<<"Start of Debugging in "<<__func__<<std::endl;
    if( (std::abs(CP_sym)<1e-12) and (std::abs(CP_brk)<1e-12))
    {
        res[0] = 0 ;
        res[1] = 0 ; 
        res[2] = 0 ; 
        return res;
    }
    double delta_theta = CP_brk-CP_sym;
    //Calculation of the thermal tanbeta
//    double omega_1_brk = gen_fluid::vcritical.at(1);
//    double omega_2_brk  = gen_fluid::vcritical.at(2);
//    double omega_CP_brk = gen_fluid::vcritical.at(3);
//    double tanbetasqrt_bar = (std::pow(omega_CP_brk,2) + std::pow(omega_2_brk,2))/(std::pow(omega_1_brk,2));
//TODO: Check for the Tbeta suppression? Does this need to be taken into account?!
    // double theta = (CP_brk - delta_theta/2. * ( 1 + std::tanh(z/gen_fluid::LW)))/(1+tanbetasqrt_bar);
    double theta = (CP_brk - delta_theta/2. * ( 1 + std::tanh(z/gen_fluid::LW)));
    // double theta_prime = - (delta_theta)/(2*gen_fluid::LW*std::pow(std::cosh(z/gen_fluid::LW),2)*(1+tanbetasqrt_bar));
    double theta_prime = - (delta_theta)/(2*gen_fluid::LW*std::pow(std::cosh(z/gen_fluid::LW),2));
    // double theta_2prime = (delta_theta*std::tanh(z/gen_fluid::LW))/(std::pow(gen_fluid::LW*std::cosh(z/gen_fluid::LW),2)*(1+tanbetasqrt_bar));
    double theta_2prime = (delta_theta*std::tanh(z/gen_fluid::LW))/(std::pow(gen_fluid::LW*std::cosh(z/gen_fluid::LW),2));

    res[0] = theta;
    res[1] = theta_prime;
    res[2] = theta_2prime;
    if(debug)
    {
//    std::cout<<sep << "tanbetasqrt_bar = " <<tanbetasqrt_bar << std::endl;
    std::cout<< sep << "Sym phase = " << CP_sym<<std::endl;
    std::cout<< sep << "Brk phase = " <<CP_brk<<std::endl;
    std::cout<< sep << "delttheta = " <<delta_theta<<std::endl;
    std::cout<< sep << "theta(z)  = " <<res[0]<<std::endl;
    std::cout<< sep << "theta_prime= " << res[1]<<std::endl;
    std::cout<< sep << "theta_2prime = " << res[2]<<std::endl;
    }
    if(debug) std::cout<<"END of Debugging in " <<__func__<<std::endl;
    return res;
}


void gen_fluid::tau_func(double z, std::vector<double>& m_lep, std::vector<double>& m_lep_prime)
{
    /**
     * Function: Calculation of the tau mass and its derivatives in respect of the VEV's
     * Return: Pass by Ref in m_lep and m_lep_prime
     */
    bool debug = false;
    if(debug) std::cout<<"Begin debug in function "<<__func__<<std::endl;
    int nlep = gen_fluid::modelPointer->get_NLepton();
//Configure z-dep VEV
    std::vector<double>gen_phi;
    std::size_t gen_nvev = gen_fluid::gen_vcritical.size();
    for(size_t i =0;i<gen_nvev;i++) gen_phi.push_back(gen_fluid::gen_vcritical.at(i)/2*(1-std::tanh(z/gen_fluid::LW)));
//Calculation of the mass mtau
    std::vector<double> mi;
    std::vector<double> grad_mtau;
    grad_mtau.clear();

    for(size_t i=1;i<=gen_nvev;i++)
    {
        mi.clear();
        mi=gen_fluid::modelPointer->LeptonMassesSquared(gen_phi , i);
        /*
            m_i: nlep eigenvalues and nlep derivatives with respect to field i
            Since the tau lepton is the heaviest lepton
            --> mtau = m_i.at(N)
        */
       grad_mtau.push_back(mi.at(2*nlep-1));
    }
    double mtau = std::sqrt(std::abs(mi.at(nlep-1)));
    double res =0;
    std::vector<double> grad_omega = gen_fluid::omegaprime(z);
    for(size_t i=0;i<grad_mtau.size();i++) res += grad_mtau.at(i)*grad_omega.at(i);
    if(std::isnan(res)) throw std::runtime_error("NaN @ tau_func");
    m_lep.resize(1);
    m_lep_prime.resize(1);
    m_lep[0] = mtau;
    m_lep_prime[0] = res;
}



void gen_fluid::set_class
    (
    const int bottom_mass_inp,
    struct GSL_integration_mubl& container,
    const Calc_Gam_M& Calc_Gam_inp,
    const Calc_Scp& Calc_Scp_inp,
    const Calc_kappa_t& Calc_kappa_inp
    )
{

    //   /*
	// function: Setting Class parameters for the ODE coeffecients
    // */
    bool debug = false;
    if(debug) std::cout<<"Begin of debug in "<<__func__<<std::endl;
    gen_fluid::bot_mass_flag    = bottom_mass_inp;
    gen_fluid::modelPointer     = container.getModelPointer();
    gen_fluid::par              = modelPointer->get_parStored();
    gen_fluid::symmetric_CP_violating_phase = container.getSymmetricCPViolatingPhase();
    gen_fluid::broken_CP_violating_phase = container.getBrokenCPViolatingPhase();
    gen_fluid::TOP_symmetric_CP_violating_phase = container.getSymmetricCPViolatingPhase_top();
    gen_fluid::BOT_symmetric_CP_violating_phase = container.getSymmetricCPViolatingPhase_bot();
    gen_fluid::TOP_broken_CP_violating_phase = container.getBrokenCPViolatingPhase_top();
    gen_fluid::BOT_broken_CP_violating_phase = container.getBrokenCPViolatingPhase_bot();
    gen_fluid::TAU_broken_CP_violating_phase = container.getBrokenCPViolatingPhase_tau();
    gen_fluid::TAU_symmetric_CP_violating_phase = container.getSymmetricCPViolatingPhase_tau();

    gen_fluid::vcritical        = container.getVEVCritical();
    gen_fluid::Calc_Gam_obj     = Calc_Gam_inp;
    gen_fluid::Calc_Scp_obj     = Calc_Scp_inp;
    gen_fluid::Calc_kappa_obj   = Calc_kappa_inp;
    gen_fluid::tanbeta          = par.at(7);
    gen_fluid::Yuk_Type         = par.at(8);
    gen_fluid::gen_vcritical=gen_fluid::modelPointer->MinimizeOrderVEV(container.getVEVCritical());
    gen_fluid::vsym     = container.get_vev_sym_theta();
    gen_fluid::LW       = container.getLW();
    gen_fluid::Temp     = container.getTC();
    gen_fluid::vw       = container.getvw();
    //Consistent with hep-ph/9410281
    gen_fluid::Dq       = 6./Temp;
    gen_fluid::Dt       = 6./Temp;
    gen_fluid::Dh       = 100./Temp;
    gen_fluid::Dlep     = 380./Temp;
    gen_fluid::Dtau     = 100./Temp;
    double sinbeta      = std::sin(std::atan(tanbeta));
    double cosbeta      = std::cos(std::atan(tanbeta));
    gen_fluid::yuk_q    = std::sqrt(2)*mtop_0/(C_vev0*sinbeta);
    gen_fluid::Gam_SS   = 14* std::pow(alphaS,4)*Temp;

    double yuk_t=0,yuk_b=0,yuk_tau=0;
    //top thermal mass
    yuk_t        = std::sqrt(2)*mtop_0/(C_vev0*sinbeta);
    auto top_thermal = Calc_ThermalMass_q(yuk_t,Temp);
    msqrt_thermal_top = top_thermal.first;
    dmsqrt_thermal_top = top_thermal.second;

    if(Yuk_Type == 1)
    {
    //bot&tau thermal mass
      yuk_b = std::sqrt(2)*mbot_0/(C_vev0*sinbeta);
      auto bot_thermal = Calc_ThermalMass_q(yuk_b,Temp);
      msqrt_thermal_bot=bot_thermal.first;
      dmsqrt_thermal_bot=bot_thermal.second;

      yuk_tau = std::sqrt(2)*mtau_0/(C_vev0*sinbeta);
      auto tau_thermal = Calc_ThermalMass_l(yuk_tau,Temp);
      msqrt_thermal_tau=tau_thermal.first;
      dmsqrt_thermal_tau=tau_thermal.second;
    } 
    if(Yuk_Type == 2)
    {
        yuk_b = std::sqrt(2)*mbot_0/(C_vev0*cosbeta);
        auto bot_thermal = Calc_ThermalMass_q(yuk_b,Temp);
        msqrt_thermal_bot=bot_thermal.first;
        dmsqrt_thermal_bot=bot_thermal.second;

        yuk_tau = std::sqrt(2)*mtau_0/(C_vev0*cosbeta);
        auto tau_thermal = Calc_ThermalMass_l(yuk_tau,Temp);
        msqrt_thermal_tau=tau_thermal.first;
        dmsqrt_thermal_tau=tau_thermal.second;
    }
    if((Yuk_Type!=1) and (Yuk_Type!=2))
    {
        throw std::runtime_error("No valid Yuk Type in gen_fluid::set_class\n"); 
    }
    
    //ArxiV: hep-ph/9410281 @ Appendix C+D WITHOUT HIGGS EXCHANGE! (only Higgs induced rates!)
    gen_fluid::Gam_Y    = 0.19*alphaS*yuk_q*yuk_q*Temp;
    gen_fluid::Gam_Y_t  = 0.19*alphaS*yuk_t*yuk_t*Temp;
    gen_fluid::Gam_Y_b  = 0.19*alphaS*yuk_b*yuk_b*Temp;
    gen_fluid::Gam_Y_tau= 0.28*alphaW*yuk_tau*yuk_tau*Temp;
    // gen_fluid::Gam_Y_u  = 0.19*alphaS*yuk_u*yuk_u*Temp;
    

    if(debug)
    {
        std::cout<<sep << "LW = " << LW<<std::endl;
        std::cout<<sep << "Temp= "<<Temp<<std::endl;
        std::cout<<sep << "tanbeta= "<<tanbeta<<std::endl;
        std::cout<<sep << "GamY_t = " <<Gam_Y_t<<std::endl;
        std::cout<<sep << "GamY_b = " <<Gam_Y_b<<std::endl;
        std::cout<<sep << "GamY_tau = " <<Gam_Y_tau<<std::endl;
        std::cout<<sep << "vsym:"<<std::endl;
        for(auto x : vsym) std::cout<<sep << sep<<x<<std::endl;
        std::cout<<sep << "vcritical:"<<std::endl;
        for(auto x : vcritical) std::cout<<sep << sep<<x<<std::endl;
    }
}

std::pair<double,double> gen_fluid::Calc_ThermalMass_q(double & YukCoupling,double & T_in){
    double msqrt = std::pow(gprime,2)/18 + std::pow(gS,2)/6 + std::pow(YukCoupling,2)/8;
    msqrt = msqrt*std::pow(T_in,2);
    double dmsqrt = 5*std::pow(gprime,2)/96 - 3*std::pow(g,2)/32 + std::pow(YukCoupling,2)/16;
    dmsqrt= dmsqrt*std::pow(T_in,2);
    return std::make_pair(msqrt,dmsqrt);
}
std::pair<double,double> gen_fluid::Calc_ThermalMass_l(double & YukCoupling,double & T_in){
    double msqrt = std::pow(gprime,2)/8 + std::pow(YukCoupling,2)/8;
    msqrt = msqrt*std::pow(T_in,2);
    double dmsqrt = 3*std::pow(gprime,2)/4 - 3*std::pow(g,2)/32 + std::pow(YukCoupling,2)/16;
    dmsqrt= dmsqrt*std::pow(T_in,2);
    return std::make_pair(msqrt,dmsqrt);
}



}//END:Namespace:Baryo
}//END:Namespace:BSMPT
