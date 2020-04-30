#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/minimizer/Minimizer.h>
#include <BSMPT/minimizer/MinimizePlane.h>

#include <BSMPT/Kfactors/Kfactors.h>
#include <BSMPT/baryo_calculation/transport_equations.h>

////////Included Models for the transport equations////////
#include <BSMPT/baryo_calculation/Fluid_Type/gen_calc.h>
#include <BSMPT/baryo_calculation/Fluid_Type/gen_func_fluid.h>
#include <BSMPT/baryo_calculation/Fluid_Type/top_source.h>
#include <BSMPT/baryo_calculation/Fluid_Type/bot_source.h>
#include <BSMPT/baryo_calculation/Fluid_Type/tau_source.h>
//////////////////////////////////////////////////////////


using namespace boost::numeric::odeint;

///////////////OPTIONS//////////////////
	/*
		method_transport == 0 --> Every method is calculated
		method_transport == 1 --> top only included in transport equations
		method_transport == 2 --> top+bot included in transport equations
		method_transport == 3 --> top+bot+tau included in transport equations
	*/
	const int method_transport = 0 ;

	/*
    	bot_mass_flag ==1 --> non-zero bot mass
    	bot_mass_flag ==2 --> zero bot mass 
	*/
	const int bot_mass_flag = 1 ;
///////////////////////////////////////

int main(int argc, char *argv[]) try{


  bool debug = true;

	int Model=-1;

	if(!( argc == 7 or argc == 8) )
	{
		std::cerr << "./Calc_eta_fluid Model Inputfile Outputfile  LineStart LineEnd v_W  \n";
		ShowInputError();
		return EXIT_FAILURE;
	}


    Model=ModelID::getModel(argv[1]);
	if(Model==-1) {
		std::cerr << "Your Model parameter does not match with the implemented Models." << std::endl;
		ShowInputError();
		return EXIT_FAILURE;
	}

	double LineStart,LineEnd;
	char* in_file;char* out_file;

	in_file = argv[2];
	out_file = argv[3];

	if(Model == -1) Model = atoi(argv[1]);
	LineStart = atoi(argv[4]);
	LineEnd = atoi(argv[5]);
  	double vw = atof(argv[6]);

	bool TerminalOutput = false;
	if(argc == 8) {
		std::string s7 = argv[7];
		std::cout << s7 << std::endl;
		TerminalOutput = ("y" == s7);
	}

	if(LineStart < 1)
	{
		std::cout << "Start line counting with 1" << std::endl;
		return EXIT_FAILURE;
	}
	if(LineStart > LineEnd)
	{
		std::cout << "LineEnd is smaller then LineStart " << std::endl;
		return EXIT_FAILURE;
	}

	int linecounter = 1;
	std::ifstream infile(in_file);
	if(!infile.good()) {
		std::cout << "Input file not found " << std::endl;
		return EXIT_FAILURE;
	}

	std::ofstream outfile(out_file);
	if(!outfile.good())
	{
		std::cout << "Can not create file " << out_file << std::endl;
		return EXIT_FAILURE;
	}
	std::string linestr;
    std::shared_ptr<Class_Potential_Origin> modelPointer = ModelID::FChoose(Model);
	// int Type;
	// double tmp;
	int nPar,nParCT;
	nPar = modelPointer->nPar;
	nParCT = modelPointer->nParCT;
	size_t ndim = modelPointer->nVEV;
	std::vector<double> par(nPar);
	std::vector<double> parCT(nParCT);
	//Declare the vector for the PTFinder algorithm
    std::vector<double> sol;

	//Initialise the polarisation for the KFACTORs ---> Not needed in the fluid type ansatz
	Kfactors_interpolated_GSL Kfunctions;
	Kfunctions.init();

	/*
		Declaration of all needed parameters
	*/
    double LW =0;


	if(debug) std::cout << "Begin of Inputread" <<std::endl;
	while(getline(infile,linestr))
	{
		if(linecounter > LineEnd) break;

		if(linecounter == 1)
		  {
		    // outfile << linestr << sep << modelPointer->addLegendCT() << sep;
		    // outfile << modelPointer->addLegendTemp();
			// outfile << std::endl;

			//Assuming we already read a BSMPT Output file 
			modelPointer->setUseIndexCol(linestr);
			outfile << linestr;
			outfile << sep <<"LW_J"<<sep<<"vw_J"<<sep<<"eta_J";
			outfile << std::endl;

		  }
		if(linecounter >= LineStart and linecounter <= LineEnd and linecounter != 1)
		{
			if(TerminalOutput)
			{
				std::cout << "Currently at line " << linecounter << std::endl;
			}
			std::pair<std::vector<double>,std::vector<double>> parameters = modelPointer->initModel(linestr);
			par=parameters.first;
			parCT = parameters.second;
			modelPointer->FindSignSymmetries();

			if(LineStart == LineEnd ) modelPointer->write();
			sol.clear();

      if(debug) std::cout<< "End of parameter inputread -> PTFinder_gen_all called"<<std::endl;
		
			PTFinder_gen_all(modelPointer,0,300,sol,3);
		//Define parameters for eta and vw
			double eta = 0;
		// Struct Container for the parameters needed for the calculation
			struct GSL_integration_mubl GSL_integration_mubl_container;

      if(TerminalOutput) std::cout << "PFTinder_gen_all END" << std::endl;
//START: SFOEWPT is found at given line 
		if(sol.at(2) == 1 and C_PT*sol.at(0) < sol.at(1) )
		{
	//Calculation of LW 
			std::vector<double> vcritical,vbarrier;
			for(size_t i=3;i<ndim+3;i++) vcritical.push_back(sol.at(i));
			int CPoddIndex=3, CPevenIndex=2;
			double brokenPhase = std::atan2(vcritical.at(CPoddIndex),vcritical.at(CPevenIndex));
		//Sign Check for the broken phase
			if(std::abs(brokenPhase) > 0.5*M_PI)
			{
				for(auto symmetry : modelPointer->SignSymmetries)
				{
					double tmpphase=std::atan2(symmetry.at(CPoddIndex)*vcritical.at(CPoddIndex), symmetry.at(CPevenIndex)*vcritical.at(CPevenIndex) );
					if(std::abs(tmpphase) <= 0.5*M_PI)
					{
						for(size_t i=0;i<vcritical.size();i++) vcritical.at(i) *= symmetry.at(i);
						break;
					}
				}
			}
			double TC = sol.at(0);
		//Symmetric Phase (symmetric minimum)-->Calculated at TC+1 
			std::vector<double> vevsymmetricSolutiontmp,checksym, startpoint;
			for(size_t i=0;i<modelPointer->nVEV;i++) startpoint.push_back(0.5*vcritical.at(i));
			Minimize_gen_all(modelPointer,TC+1,vevsymmetricSolutiontmp,checksym,startpoint);

			if(debug)
			{
				std::cout << "Solution for the symmetric phase  : " << std::endl;
				for(size_t i=0;i<vevsymmetricSolutiontmp.size();i++)
				{
					std::cout << "vevsymmetricSolution[" << i << "] = " << vevsymmetricSolutiontmp.at(i) << std::endl;
				}
			}
			double absVEVSymmetricSolutionTmp = 0;
			for(auto x:vevsymmetricSolutiontmp) absVEVSymmetricSolutionTmp+=std::pow(x,2);
			if(absVEVSymmetricSolutionTmp != 0) std::cout << "Check for sign of non vanishing components of symmetric vacuum : " << std::endl;
//TODO Implementation for the CN2HDM ?!
		//Initialisation for the parameter container
			if(debug) std::cout<<"Start init container"<<std::endl;
			//init needs vevsymmetricSolutiontmp! The real vevconfig for the actual top phase is calculated in init() --> vev_sym_theta
			GSL_integration_mubl_container.init(vw,vcritical,vevsymmetricSolutiontmp,TC,modelPointer,Kfunctions);
			if(debug) std::cout<<"End init container"<<std::endl;
			GSL_integration_mubl_container.setZMAX(5,true);
			if(debug) std::cout<<"zmax = " << GSL_integration_mubl_container.getZMAX()<<std::endl;
			

		//START: mu_bl calculation
            int n_step=50;
            std::vector<double> vec_mubl;
            std::vector<double> vec_z;
            vec_mubl.resize(n_step);
            vec_z.resize(n_step);

            if(TerminalOutput) std::cout<<"___Set Up nL Grid___ "<<std::endl;
		//Initialisation of Gam_M class for the numerical evaluation of Gamma_M
			Calc_Gam_M Calc_Gam_inp;
		//Initialisation of S_CP class for the numerical evaluation of S_CP
            Calc_Scp Calc_Scp_inp;
		//Initialisation of kappa class for the numerical evaluation of kappa
            Calc_kappa_t Calc_kappa_inp;

		//Initialisation of fluid transport class for the numerical evaluation of the transport equations
		
		//TOP only
            top_source C_top;
		    C_top.set_class(bot_mass_flag, GSL_integration_mubl_container,Calc_Gam_inp,Calc_Scp_inp,Calc_kappa_inp);
		// //BOT+TOP	
			bot_source C_bot;
		    C_bot.set_class(bot_mass_flag, GSL_integration_mubl_container,Calc_Gam_inp,Calc_Scp_inp,Calc_kappa_inp);
		// //TAU+BOT+TOP
			tau_source C_tau;
		    C_tau.set_class(bot_mass_flag, GSL_integration_mubl_container,Calc_Gam_inp,Calc_Scp_inp,Calc_kappa_inp);
		//Decide which method 
			std::pair<std::vector<double> , std::vector<double> > arr_nL;
			switch(GSL_integration_mubl_container.get_transport_method())
			{
				case 1:  arr_nL = set_up_nL_grid(n_step , GSL_integration_mubl_container , C_top );
						break;
				case 2:  arr_nL = set_up_nL_grid(n_step , GSL_integration_mubl_container , C_bot );
						break;
				case 3:  arr_nL = set_up_nL_grid(n_step , GSL_integration_mubl_container , C_tau );
						break;
				default: throw std::runtime_error("No valid transport method chosen \n");
				
			}
			
		// 	//Initialisation of the eta class for the numercial evaluation of eta 
            Calc_eta C_eta;
            C_eta.set_class(arr_nL,TC,vw);
			eta = Nintegrate_eta(C_eta,0 , GSL_integration_mubl_container.getZMAX());
			if(TerminalOutput) std::cout<<"eta = n_b/s = "<<eta<<std::endl;
          }//END: SFOEWPT
			/*
				Start of outfile 
			*/
		outfile << linestr;
		outfile << sep<<LW<<sep<<vw<<sep<<eta<<std::endl;
       }//END: Valid Line Number
		linecounter++;
		if(infile.eof()) break;
	}//END: Line Calculation
	if(TerminalOutput) std::cout << std::endl;




	outfile.close();


  return EXIT_SUCCESS;

}
catch(std::exception& e){
		std::cerr << e.what() << std::endl;
		return EXIT_FAILURE;
}
