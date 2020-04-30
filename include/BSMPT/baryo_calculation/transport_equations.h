/*
 * transport_equations.h
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
 #ifndef SRC_transport_equations_H_
 #define SRC_transport_equations_H_


/**
 * @file
 */


 #include <BSMPT/models/IncludeAllModels.h>
 #include <boost/numeric/odeint.hpp>
 #include <boost/numeric/odeint/iterator/const_step_time_iterator.hpp>
 #include <boost/math/interpolators/cubic_b_spline.hpp>




namespace BSMPT{
  class Class_Potential_Origin;
//  namespace Kfactors {
//    class Kfactors_interpolated_GSL;
//  }

  namespace Baryo{



  typedef std::vector< double > state_type;
  typedef boost::numeric::odeint::runge_kutta_cash_karp54< state_type > error_stepper_type;



  /**
  * This class handles the evaluation of the transport equations as shown in arXiv:hep-ph/0605242v2 Eq (44) and (45)
  * with the thermal velocity. The transport equations with second order ODEs in the chemical potentials is given in
  * TODO:: paper/thesis
  */
  class transport_equations {
  protected:
      /**
       * @brief UseTanBetaSuppression Use the thermal tanbeta suppression in the calculation of the theta(z)
       */
      bool UseTanBetaSuppression = false;
      /**
      * If true the transport equations with the plasma velocitys are used, otherwise the second order transport equations
      * in the chemical potentials are used.
      */
      bool UseVelocityTransportEquations=false;
      /**
      * Internal storage for the wall velocity
      */
      double vw;
      /**
      * Internal storage for the wall thickness
      */
      double LW;
      /**
      * Internal storage for the critical temperature
      */
      double TC;
      /**
      * Internal storage for the CP violating phase of the symmetric minimum
      */
      double symmetric_CP_violating_phase = -500000;
      /**
      * CP violating phase in the broken vacuum
        */
      double broken_CP_violating_phase = -500000;
      /**
      * Pointer to the specific model under investigation. At this moment the transport equations are only implemented and
      * tested for the C2HDM.
      */
      std::shared_ptr<Class_Potential_Origin> modelPointer;

      /**
      * Internal storage for the value of the VEVs at the broken minimum at the critical temperature.
      */
      std::vector<double> vev_critical;



  public:
      /**
       * @brief transport_equations
       * @param params A GSL_integration_mubl struct with the information used
       */
      transport_equations(const struct GSL_integration_mubl& params);
      virtual ~transport_equations();


      /**
    * Overloads the () operator to calculate the ODE for the transport equation
    * @param x Depending on transport_equations::UseVelocityTransportEquations the value of the chemical potentials and either their
    * derivatives or the thermal velocities at the point z.
    * @param dxdt Here the derivative of x at the point z will be stored
    * @param z The distance to the wall (z = + infinity is the symmetric minimum, z = - infinity is the broken minimum and
    * z=0 is the wall)
    */
      transport_equations &operator()	( const state_type &x , state_type &dxdt , const double /* z */ );

      /**
    * Calculates the top mass and its derivative w.r.t the single components of the VEVs.
    * @param vev The VEV configuration at which the mass and the derivatives should be calculated.
    *
    * @return vector which will store the result. Will be cleared in the process of the function call. This will have 9 entries, the first
    is the mass and then the 8 derivatives.
    */
      std::vector<double> get_top_mass_and_derivative(const std::vector<double>& vev) const;

      /**
    * Calculates the W mass at a given VEV and temerpature. As we assume no charge breaking vev in the c2hdm we look for the
    * degenerate eigenvalue in the gauge boson mass matrix.
    * @param vev The VEV configuration at which to compute the W mass
    * @param T The temperature at which to compute the W mass
    */
      double get_W_mass(const std::vector<double>& vev, const double& T) const;

      /**
    * Calculates the VEV at a given distance z from the wall
    * @param z Distance to the wall. z < 0 is inside the broken phase and z > 0 is in the symmetric phase.
    * @return vector in which the VEV configuration will be stored
    */
      std::vector<double> calculate_vev(const double& z) const;

      /**
    * Calculates the derivatives of the VEV at a given distance z from the wall
    * @param z Distance to the wall. z < 0 is inside the broken phase and z > 0 is in the symmetric phase.
    * @return vector in which the VEV configuration will be stored
    */
      std::vector<double> calculate_vev_derivative(const double& z) const;

      /**
    * Calculates the CP violating angle theta according to 0605242 and its derivatives
    * @param z Distance to the wall. z < 0 is inside the broken phase and z > 0 is in the symmetric phase.
    * @param diff Switch if the actual value should be returned (0), the first derivative (1) or the second (2)
    */
      double calculate_theta(const double& z, const int& diff) const;


  };


  /**
 * Struct for the GSL integration over mu_BL
 */
  struct GSL_integration_mubl{
  private:
      /**
          * If true the transport equations with the plasma velocitys are used, otherwise the second order transport equations
          * in the chemical potentials are used.
          */
      bool UseVelocityTransportEquations=false;

      /**
        * model information
        */
      std::shared_ptr<Class_Potential_Origin> modelPointer;
      /**
       * wall velocity
       */
      double vw;

      /**
        * wall thickness
        */
      double LW;

      /**
        * maximal value at which mu = 0 is assumed
        */
      double zmax;

      /**
        * temperature
        */
      double TC;

      /**
        * EW VEV at TC
        */
      double vc;

      /**
        * vev in the critical phase
        */
      std::vector<double> vev_critical;
      /**
        * vev in the symmetric phase
        */
      std::vector<double> vev_symmetric;

      /**
        * CP violating phase in the symmetric vacuum
        */
      double symmetric_CP_violating_phase = -500000;

      /**
        * CP violating phase in the broken vacuum
        */
      double broken_CP_violating_phase = -500000;

      /**
         * CP violating phase of the top quark in the symmetric vacuum
         */
      double TOP_symmetric_CP_violating_phase;
      /**
         * CP violating phase of the top quark in the broken vacuum
         */
      double TOP_broken_CP_violating_phase;
      /**
         * CP violating phase of the bot quark in the symmetric vacuum
         */
      double BOT_symmetric_CP_violating_phase;
      /**
         * CP violating phase of the bot quark in the broken vacuum
         */
      double BOT_broken_CP_violating_phase;
      //TODO: Calculation of the Tau phase
      /**
        * CP violating phase of the tau lepton in the symmetric vacuum
        */
      double TAU_symmetric_CP_violating_phase;
      /**
        * CP violating phase of the tau quark in the broken vacuum
        */
      double TAU_broken_CP_violating_phase;
      /**
         * Vector with parameter of the given point
         */
      std::vector<double> par;
      /**
         * Method for the transport equations
         * 1 --> top
         * 2 --> bot
         * 3 --> tau
         */
      int transport_method =  1 ;


      /**
         * VEV configuration before symmetric phase (for the top phase calculation)
         */
      std::vector<double> vev_sym_theta;


  public:
      /**
        * set the CP violating phase in the symmetric vacuum
        */
      void setSymmetricCPViolatingPhase(double Phase);
      /**
        * get the CP violating phase in the symmetric vacuum
        */
      double getSymmetricCPViolatingPhase() const;
      /**
        * get the CP violating phase in the symmetric vacuum
        */
      double getBrokenCPViolatingPhase() const;
      /**
        * get the CP-violating phase of the top quark in the symmetric vacuum
        */
      double getSymmetricCPViolatingPhase_top() const;
      /**
        * get the CP-violating phase of the bot quark in the symmetric vacuum
        */
      double getSymmetricCPViolatingPhase_bot() const;
      /**
        * get the CP-violating phase of the tau quark in the symmetric vacuum
        */
      double getSymmetricCPViolatingPhase_tau() const;

      /**
        * get the CP-violating phase of the top quark in the broken vacuum
        */
      double getBrokenCPViolatingPhase_top() const;
      /**
        * get the CP-violating phase of the bot quark in the broken vacuum
        */
      double getBrokenCPViolatingPhase_bot() const;
      /**
        * get the CP-violating phase of the tau quark in the broken vacuum
        */
      double getBrokenCPViolatingPhase_tau() const;




      /**
        * get critical VEV
        */
      std::vector<double> getVEVCritical() const;
      /**
        * get vev_symmetric
        */
      std::vector<double> getVEVsym() const;
      /**
        * get TC
        */
      double getTC() const;
      /**
        * set TC
        */
      void setTC(double TC_in);
      /**
        * get zmax
        */
      double getZMAX() const;
      /**
        * set vw
        */
      void setvw(double vw_in);
      /**
        * get vw
        */
      double getvw() const;
      /**
        * @brief getLW
        * @return LW
        */
      double getLW() const;

      /**
       * @brief setpar sets par
       */
      void setpar(std::vector<double>);
      /**
       * @brief getpar
       * @return par
       */
      std::vector<double> getpar() ;
      /**
         * Set function to chose the method for the transport equations
        */
      void set_transport_method(int method);
      /**
       * @brief get_transport_method
       * @return transport_method
       */
      int get_transport_method();
      /**
       * @brief setZMAX defines the value to treat mu(ZMAX) = 0
       * @param z_in new value to set zMAX to
       * @param MultiplesOfLW if true sets zMAX = z_in * LW
       */
      void setZMAX(double z_in,bool MultiplesOfLW)  ;

      /**
       * @brief set_vev_sym_theta
       * @param vev_in new value of vev_sym_theta
       */
      void set_vev_sym_theta(std::vector<double> & vev_in);
      /**
       * @brief get_vev_sym_theta
       * @return  vev_sym_theta
       */
      std::vector<double> get_vev_sym_theta() const ;


      /**
        * set the UseVelocityTransportEquations parameter
        */
      void setUseVelocityTransportEquations(bool in);

      /**
        * @brief getUseVelocityTransportEquations
        * @return the UseVelocityTransportEquations parameter
        */
      bool getUseVelocityTransportEquations() const;


      /**
       * @brief init initialises the parameters of struct
       * @param vw_input new value of vw
       * @param vev_critical_input new value of vev_critical
       * @param vev_symmetric_input new value of vev_symmetric
       * @param TC_input new value of TC
       * @param modelPointer_input shared_ptr for the model
       */
      void init(const double& vw_input, std::vector<double>& vev_critical_input, std::vector<double>& vev_symmetric_input,
                const double& TC_input,
                std::shared_ptr<Class_Potential_Origin>& modelPointer_input);

      /**
       * @brief getModelPointer
       * @return the shared_ptr to the model
       */
      std::shared_ptr<Class_Potential_Origin> getModelPointer() const;
  };

  /**
  * struct for the integration of the 1D interpolation of mu_BL
  */
  struct GSL_mubl_interpolation{
      /**
      * Boost interpolation of mu_BL, calculated through generate_mubl_spline
      */
      boost::math::cubic_b_spline<double> spline;
      /**
      * wall velocity
      */
      double vw;
      /**
      * Temperature
      */
      double TC;
  };

  /**
  * Solve the Transport equations using the transport_equations class.
  * @param z Distance from the wall (z=0) at which the equations should be solved. The symmetric phase is at z=infinity.
  * @param parStart The boundary conditions in the symmetric phase for the chemical potentials
  * @param params The GSL_integration_mubl struct which contains all necessary informations of the model.
  * @return The values of the chemical potentials at distance z will be stored
  */
  std::vector<double> calculateTransportEquation(
          const double& z,
          const std::vector<double>& parStart,
          const struct GSL_integration_mubl& params);

  /**
  * Evaluates 0605242 Eq (47) at distance z from the wall
  * @param z distance from the wall
  * @param p GSL_integration_mubl object with the necessary information for the integration
  */
  double mubl_func(double z, void *p);


  /**
  * Calculates the integrand necessary for the baryon-antibaryon asymmetry
  * @param z distance to the wall
  * @param p void pointer to a GSL_integration_mubl struct
  */
  double eta_integrand_func(double z,  void *p);

  /**
  * Integrate over mu_BL without using the interpolation for mu_BL
  * @param p GSL_integration_mubl without the mubl spline
  */
  double Integrate_mubl(const struct GSL_integration_mubl& p);

  /**
  * generate the spline of 0605242 Eq (47)
  * @param p GSL_integration_mubl object which containts the necessary information about the integration
  * @param nstep number of steps between the wall (z=0) and the symmetric phase (z=3 LW )
  * @return GSL_mubl_interpolation object in which the created spline will be stored
  */
  GSL_mubl_interpolation generate_mubl_spline(const struct GSL_integration_mubl& p, int nstep);

  /**
  * evaluates the spline and multiples it with the prefacotr exp(-nu z) in 0605242 Eq (48)
  * @param z distance from the wall
  * @param p GSL_mubl_interpolation object with the spline of mu_BL calculated in generate_mubl_spline
  */
  double mubl_interpolation(double z,void *p);

  /**
  * Integrates over mubl_interpolation and evaluates 0605242 Eq (48)
  * @param p GSL_integration_mubl with the mubl spline
  */
  double Integrate_mubl_interpolated(const struct GSL_integration_mubl& p);




  /**
 * struct for intermediate outputs during the ODE solving.
 */
  struct push_back_state_and_time
  {
      std::vector< state_type >& m_states;
      std::vector< double >& m_times;

      push_back_state_and_time( std::vector< state_type > &states , std::vector< double > &times )
          : m_states( states ) , m_times( times ) { }

      void operator()( const state_type &x , double t )
      {
          m_states.push_back( x );
          m_times.push_back( t );
      }
  };


  }
}

#endif /* SRC_transport_equations_H_ */

