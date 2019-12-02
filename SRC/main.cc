#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include "ConfigFile.tpp"
using namespace std;

double randomDouble(const double& min, const double& max); // function ton generate a random (uniform) double between 'min' and 'max'
double potentialAction(const double& m, const double& omega, const vector<double>& path, const double& displacement=0.0);
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

  string inputPath("configuration.in"); // Fichier d'input par defaut
  if(argc>1) // Fichier d'input specifie par l'utilisateur ("./Exercice7 config_perso.in")
    inputPath = argv[1];

  ConfigFile configFile(inputPath); // Les parametres sont lus et stockes dans une "map" de strings.

  for(int i(2); i<argc; ++i) // Input complementaires ("./Exercice7 config_perso.in input_scan=[valeur]")
    configFile.process(argv[i]);

  unsigned int N_part(configFile.get<unsigned int>("N_part")), // number of particles in the system
               N_slices(configFile.get<unsigned int>("N_slices")); // number of (imaginary) time slices
  double tau_final(configFile.get<double>("tau")); // tau<->beta=1/(T*k_B) the so called inverse temperature
	double d_tau(tau_final/N_slices); // imaginary time step
	double m(configFile.get<double>("mass")*d_tau), // dimensionless mass
         omega(configFile.get<double>("frequency")*d_tau); // dimensionless angular frequency

  string type_V(configFile.get<string>("type_V"));
	Potential* V;
	if(type_V=="harmonic") V = new Potential_harm(configFile);
	else if(type_V=="double") V = new Potential_double(configFile);
	else if(type_V=="square") V = new Potential_square(configFile);
  else{
    cerr << "Please choose between ""harmonic"", ""double"" or ""square"" for ""type_V""." << endl;
    return -1;
  }
  double pos_min(configFile.get<double>("pos_min")), // initial mininum position
         pos_max(configFile.get<double>("pos_max")); // initial maximal position
  double h(configFile.get<double>("h")); // initial maximum displacement of a point in the path
  double accrate(0.0), idrate(configFile.get<double>("idrate"));
	vector<vector<double>> system(N_part,vector<double>(N_slices,0.0)); // table : each row corresponds to a particle and contains the coordinates (1D for now) in different time slices
	vector<vector<double>> verif(N_part,vector<double>(N_slices,0.0));

	double lol_test(0.0);
	for(auto& particle : system){ // initialize random paths for each particles
		for(auto& pos : particle){
			pos=randomDouble(pos_min+lol_test,pos_max+lol_test);
		}
		//l	lol_test+=3;
	}

	unsigned int MCS(configFile.get<unsigned int>("MCS")); // number of Monte Carlo iterations
  size_t n_stride(configFile.get<size_t>("n_stride")); // output is written in the output file every n_stride iterations
	unsigned int mm(0), mm_plu(0), mm_min(0); // mm : time slice randomly selected during each iteration, mm_plu=mm+1, mm_min=mm-1;
	unsigned int nn(0); // particle randomly selected during each iteration
	double new_r(0.0), displacement(0.0); // new position proposed and displacement proposed
	double s_old(0.0), s_new(0.0); // part of the action that is changed with the new position proposed (old and new values respectively)

	string output(configFile.get<string>("output")+".out"); // output file
  ofstream fichier_output(output.c_str());
  fichier_output.precision(15);
	for(auto& particle : system){
		for(auto& pos : particle){
			fichier_output << pos << " ";
		}
	}
	fichier_output << endl;

	for(size_t i(0); i<MCS; i++){

		for(size_t j(0); j<N_part; j++){ // should we directly make one 'for i=0:N_slices*N_part' ???
			for(size_t k(0); k<N_slices; k++){
				mm=rand()%N_slices; // random integer between 0 and N_slices-1
				mm_plu=(mm+N_slices-1)%N_slices; // mm-1 with periodic boundary condition
				mm_min=(mm+1)%N_slices; // mm+1 with periodic boundary condition
				nn=rand()%N_part; // random integer between 0 and N_part-1
				verif[nn][mm]++;
				new_r=system[nn][mm]+randomDouble(-1,1)*h; // proposed new position for the particle nn at time slice mm*d_tau

				// as we take the difference of new and old action S_new-S_old, we can consider only the part of the action that is affected by the proposed new position
				s_old=0.5*m*pow((system[nn][mm_plu]-system[nn][mm])/d_tau,2)
				     +0.5*m*pow((system[nn][mm]-system[nn][mm_min])/d_tau,2)
				     //+0.5*pow(system[nn][mm],2); // harmonic potential
						 +0.1*m*omega*omega/d_tau/d_tau*pow(system[nn][mm],4)
						 -5.0*m*omega*omega/d_tau/d_tau*pow(system[nn][mm],2);
				s_new=0.5*m*pow((system[nn][mm_plu]-new_r)/d_tau,2)
				     +0.5*m*pow((new_r-system[nn][mm_min])/d_tau,2)
				     //+0.5*pow(new_r,2);
						 +0.1*m*omega*omega/d_tau/d_tau*pow(new_r,4)
						 -5.0*m*omega*omega/d_tau/d_tau*pow(new_r,2);

				if(randomDouble(0,1)<exp(s_old-s_new)){ // metropolis acceptance
					system[nn][mm]=new_r;
	        accrate+=1.0/N_slices;
				}

				nn=rand()%N_part; // random integer between 0 and N_part-1
				displacement=randomDouble(-1,1)*h; // proposed displacement of the 'entire' particle nn

				// to relative move between the time slices --> only the potential action changes
				s_old=potentialAction(m,omega/d_tau,system[nn]);
				s_new=potentialAction(m,omega/d_tau,system[nn],displacement);

				if(randomDouble(0,1)<exp(s_old-s_new)){ // metropolis acceptance
					for(auto& el : system[nn]){
						el+=displacement;
					}
	        //accrate+=1.0/N_slices;
				}
			}
			//if(i>500) h*=accrate/idrate;
			//cout << accrate << " " << h << endl;
			accrate=0.0;
		}

		// writing in the output file
		if(!(i%n_stride)){
			for(const auto& particle : system){
				for(const auto& pos : particle){
					fichier_output << pos << " ";
				}
			}
			fichier_output << endl;
		}

	}

	fichier_output.close();

	double visit_min(N_part*N_slices), visit_max(0.0);
	for(const auto& el : verif){
		for(const auto& n : el){
			//cout << n/MCS << endl;
			if(n/MCS<visit_min) visit_min=n/MCS;
			if(n/MCS>visit_max) visit_max=n/MCS;

		}
		cout << endl;
	}
	cout << "Site le moins visité, moyenne : " << visit_min << endl;
	cout << "Site le plus visité,  moyenne : " << visit_max << endl;

	return 0;
}

double randomDouble(const double& min, const double& max){
	return min + (double)rand()/RAND_MAX*(max-min);
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
