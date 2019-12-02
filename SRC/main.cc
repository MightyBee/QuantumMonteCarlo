#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <ctime>
#include "ConfigFile.tpp"
using namespace std;



/*############################## NOTES ##############################
	- theoretically x_0, x_1, ... , x_(N_slices)					, (N_slices+1) points
	- we consider boundary conditions, x_0 = x_(N_slices)
	- hence we only consider x_0, x_1, ... , x_(N_slices-1)	, N_slices points
	- to get the "full picture" simply add one more point, equal to x_0

	- units are s.t. hbar = c = 1
*/



//############################## HEADERS ##############################

// Generate a random (uniform) double between 'min' and 'max'
double randomDouble(const double& min, const double& max);

// ???
double potentialAction(const double& m, const double& omega, const vector<double>& path, const double& displacement = 0.0);

// Potential function
double V(const double& x);		// V(x)
double dV(const double& x);	// V'(x)

//double QLagrangian(const vector<vector<double>>& pos, const double& d_tau);
//double diff_QLagrangian(const vector<vector<double>>& pos, const double& d_tau, const unsigned int& m, const unsigned int& n, double const& new_pos);

class Potential {
public:
  // Methodes virtuelles pures => classe abstraite
  virtual double operator()(const double& x) const = 0; // return V at point x
};

class Potential_harm: public Potential {
public:
  // Pas de constructeur par defaut => on force a specifier une valeur
  Potential_harm(ConfigFile const& configFile) :
  Potential(), m(configFile.get<double>("mass")), omega2(pow(configFile.get<double>("omega"),2))
  {}

  // Definition des methodes virtuelles pures :
  double operator()(const double& x) const
  {
    return m*omega2*x*x;
  }

private:
  double m, omega2;
};

class Potential_double: public Potential {
public:

  Potential_double(ConfigFile const& configFile) :
  Potential(),
  V0(configFile.get<double>("V0")),
  x0(configFile.get<double>("x0"))
  {}

  double operator()(const double& x) const
  {
    return V0*pow(pow(x/x0,2)-1,2);
  }

private:
  double V0, x0;
};

class Potential_square: public Potential {
public:

  Potential_square(ConfigFile const& configFile) :
  Potential(),
  V0(configFile.get<double>("V0")),
  xc(configFile.get<double>("xc")),
  L(configFile.get<double>("L"))
  {}

  double operator()(const double& x) const
  {
		if(x<xc-0.5*L || x>xc+0.5*L){
			return 0;
		}else{
			return V0;
		}
  }

private:
  double V0, xc, L;
};



int main(int argc, char* argv[]){

	//############################## GENERATE RANDOM NUMBERS ##############################

	//Future code here. For now we just do
	srand(time(0));



	//############################## VILLARD LIBRARY ##############################

  string inputPath("configuration.in"); // Fichier d'input par defaut
  if(argc>1) // Fichier d'input specifie par l'utilisateur ("./Exercice7 config_perso.in")
    inputPath = argv[1];

  ConfigFile configFile(inputPath); // Les parametres sont lus et stockes dans une "map" de strings.

  for(int i(2); i<argc; ++i) // Input complementaires ("./Exercice7 config_perso.in input_scan=[valeur]")
    configFile.process(argv[i]);


	//############################## READ PARAMETERS ##############################

	unsigned int N_part(configFile.get<unsigned int>("N_part"));			// number of particles in the system
	unsigned int N_slices(configFile.get<unsigned int>("N_slices"));		// number of (imaginary) time slices
	unsigned int N_sweeps(configFile.get<unsigned int>("N_sweeps"));		// number of Monte Carlo iterations (aka sweeps)
	double beta(configFile.get<double>("beta"));									// inverse temperature
	double d_tau(beta/N_slices);														// imaginary time step
	double mass(configFile.get<double>("mass"));									// standard mass
	double m(mass * d_tau);																// dimensionless mass
	double omega(configFile.get<double>("frequency") * d_tau);				// dimensionless angular frequency
	string type_V(configFile.get<string>("type_V"));
	Potential* V;
	if(type_V=="harmonic") V = new Potential_harm(configFile);
	else if(type_V=="double") V = new Potential_double(configFile);
	else if(type_V=="square") V = new Potential_square(configFile);
  else{
    cerr << "Please choose between ""harmonic"", ""double"" or ""square"" for ""type_V""." << endl;
    return -1;
  }
	double pos_min(configFile.get<double>("pos_min"));							// initial minimum position
	double pos_max(configFile.get<double>("pos_max"));							// initial maximal position
	double h(configFile.get<double>("h"));											// maximum uniform displacement of a point in the path
	double accrate(0.0);																	// ???
	double idrate(configFile.get<double>("idrate"));							// ???
	size_t n_stride(configFile.get<size_t>("n_stride"));						// output is written every n_stride iterations
	// Each row corresponds to a particle and contains the coordinates in different time slices
	vector<vector<double>> system(N_part, vector<double>(N_slices, 0.0));
	vector<vector<int>> verif(N_part, vector<int>(N_slices, 0));
	//Output file
	string output(configFile.get<string>("output")+".out"); // output file
  ofstream fichier_output(output.c_str());
	fichier_output.precision(15);														// Precision



	//############################## DEFINE VARIABLES ##############################
	unsigned int mm(0), mm_plu(0), mm_min(0); // mm : time slice randomly selected during each iteration, mm_plu=mm+1, mm_min=mm-1;
	unsigned int nn(0); // particle randomly selected during each iteration
	double new_r(0.0), displacement(0.0); // new position proposed and displacement proposed
	double s_old(0.0), s_new(0.0); // part of the action that is changed with the new position proposed (old and new values respectively)



	//############################## INITIAL PATH AND OUTPUT ##############################

	for(auto& particle : system){ // initialize random paths for each particles
		for(auto& pos : particle){
			pos = randomDouble(pos_min, pos_max);
			fichier_output << pos << " ";
		}
	}
	fichier_output << endl;



	//############################## METROPOLIS ALGORITHM ##############################
	//For every sweep...
	for(size_t i(0); i < N_sweeps; i++){
		//For every particle...
		// ??? should we directly make one 'for i=0:N_slices*N_part' ???
		for(size_t j(0); j < N_part; j++){
			for(size_t k(0); k < N_slices; k++){

				// Select random particle and random point
				mm = rand()%N_slices; // random integer between 0 and N_slices-1
				mm_plu = (mm + N_slices - 1)%N_slices; // mm-1 with periodic boundary condition
				mm_min = (mm + 1)%N_slices; // mm+1 with periodic boundary condition
				nn = rand()%N_part; // random integer between 0 and N_part-1

				verif[nn][mm]++;	//this point has been visited one more time

				new_r = system[nn][mm] + h * randomDouble(-1, 1); // proposed new position

				// as we take the difference of new and old action S_new-S_old, we can
				// consider only the part of the action that is affected by the proposed new position
				s_old = 0.5*m*pow((system[nn][mm_plu]-system[nn][mm])/d_tau,2)
						+0.5*m*pow((system[nn][mm]-system[nn][mm_min])/d_tau,2)
						//+0.5*pow(system[nn][mm],2); // harmonic potential
						+0.1*m*omega*omega/d_tau/d_tau*pow(system[nn][mm],4)
						-5.0*m*omega*omega/d_tau/d_tau*pow(system[nn][mm],2);
				s_new = 0.5*m*pow((system[nn][mm_plu]-new_r)/d_tau,2)
						+0.5*m*pow((new_r-system[nn][mm_min])/d_tau,2)
						//+0.5*pow(new_r,2);
						+0.1*m*omega*omega/d_tau/d_tau*pow(new_r,4)
						-5.0*m*omega*omega/d_tau/d_tau*pow(new_r,2);

				if(randomDouble(0, 1) <= exp(-(s_new - s_old))){ // metropolis acceptance
					system[nn][mm] = new_r;		// update position with new one
					accrate += 1.0/N_slices;	// ???
				}

				// ???

				/*
				nn = rand()%N_part; // random integer between 0 and N_part-1
				displacement = h * randomDouble(-1,1); // proposed displacement of the 'entire' particle nn

				// to relative move between the time slices --> only the potential action changes
				s_old=potentialAction(m,omega/d_tau,system[nn]);
				s_new=potentialAction(m,omega/d_tau,system[nn],displacement);

				if(randomDouble(0,1)<exp(s_old-s_new)){ // metropolis acceptance
					for(auto& el : system[nn]){
						el+=displacement;
					}
					//accrate+=1.0/N_slices;
				}
				*/
			}
			//if(i>500) h*=accrate/idrate;
			//cout << accrate << " " << h << endl;
			accrate = 0.0;
		}



		//############################## OUTPUT IN FILE ##############################
		if((i%n_stride) == 0){
			for(const auto& particle : system){
				for(const auto& pos : particle){
					fichier_output << pos << " ";
				}
			}
			fichier_output << endl;
		}
	}
	fichier_output.close();



	//############################## STATS ABOUT VISITING ##############################
	//Impossible values so that we are sure to catch one point
	double visit_min(N_part * N_slices + 1);
	double visit_max(-1.0);

	for(const auto& el : verif){
		for(const auto& n : el){
			//cout << n/N_sweeps << endl;
			if(n/N_sweeps < visit_min){visit_min = n/N_sweeps;}
			if(n/N_sweeps > visit_max){visit_max = n/N_sweeps;}
		}
		cout << endl;
	}
	cout << "Least visited point, mean: " << visit_min << endl;
	cout << "Most visited point, mean: " << visit_max << endl;

	return 0;
}

//########################### FUNCTION DEFINITIONS #############################

double randomDouble(const double& min, const double& max){
	return (min + (max-min) * (double)rand()/RAND_MAX);
}

double potentialAction(const double& m, const double& omega, const vector<double>& path, const double& displacement){
	double sp(0.0);
	for(const auto& el : path){
		sp+=//0.5*pow(el+displacement,2);
			+0.1*m*omega*omega*pow(el+displacement,4)
			-5.0*m*omega*omega*pow(el+displacement,2);
	}
	return sp;
}

double V(const double& x){
	return pow(x, 2);
}

double dV(const double& x){
	return (2 * x);
}

/*
double Lagrangian(const vector<vector<double>>& pos, const double& d_tau){
	unsigned int i(0);
	return pow((pos[i+1][0]-pos[i][0])/d_tau,2)+0.5*pow(pos[i][0],2);
}
*/
/*
double diff_QLagrangian(const vector<vector<double>>& pos, const double& d_tau, const unsigned int& m, const unsigned int& n, const double& new_r){
		double Lagr_old(0.0);
		double Lagr_new(0.0);
		vector<vector<double>> pos_new(pos);
		pos_new[m][n]=new_r;
		for(size_t i(0); i<pos.size()-1; i++){
			Lagr_old+=pow((pos[i+1][0]-pos[i][0])/d_tau,2)+0.5*100*pow(pos[i][0],2);
			Lagr_new+=pow((pos_new[i+1][0]-pos_new[i][0])/d_tau,2)+0.5*100*pow(pos_new[i][0],2);
		}
		return Lagr_old-Lagr_new;
}
*/
//##############################################################################
