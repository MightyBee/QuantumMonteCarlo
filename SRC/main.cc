#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include "ConfigFile.tpp"
using namespace std;

double randomDouble(const double& min, const double& max); // function ton generate a random (uniform) double between 'min' and 'max'
//double QLagrangian(const vector<vector<double>>& pos, const double& d_tau);
//double diff_QLagrangian(const vector<vector<double>>& pos, const double& d_tau, const unsigned int& m, const unsigned int& n, double const& new_pos);

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
  double pos_min(configFile.get<double>("pos_min")), // initial mininum position
         pos_max(configFile.get<double>("pos_max")); // initial maximal position
  double h(configFile.get<double>("h")); // initial maximum displacement of a point in the path
  double accrate(0.0), idrate(configFile.get<double>("idrate"));
	vector<vector<double>> system(N_part,vector<double>(N_slices,0.0)); // table : each row corresponds to a particle and contains the coordinates (1D for now) in different time slices
	vector<vector<double>> verif(N_part,vector<double>(N_slices,0.0));

	int lol_test(0);
	for(auto& particle : system){ // initialize random paths for each particles
		for(auto& pos : particle){
			pos=randomDouble(pos_min+lol_test,pos_max+lol_test);
		}
		lol_test+=3;
	}

	unsigned int MCS(configFile.get<unsigned int>("MCS")); // number of Monte Carlo iterations
  size_t n_stride(configFile.get<size_t>("n_stride")); // output is written in the output file every n_stride iterations
	unsigned int mm(0), mm_plu(0), mm_min(0); // mm : time slice randomly selected during each iteration, mm_plu=mm+1, mm_min=mm-1;
	unsigned int nn(0); // particle randomly selected during each iteration
	double new_r(0.0); // new position proposed
	double s_old(0.0), s_new(0.0); // part of the action that is changed with the new position proposed (old and new values respectively)

	string output("output.out"); // output file
  ofstream fichier_output(output.c_str());
  fichier_output.precision(15);

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
				s_old=0.5*pow(system[nn][mm_plu]-system[nn][mm],2)
				     +0.5*pow(system[nn][mm]-system[nn][mm_min],2)
				     //+0.5*pow(system[nn][mm],2); // harmonic potential
						 +0.1*pow(system[nn][mm],4)
						 -1.0*pow(system[nn][mm],2);
				s_new=0.5*pow(system[nn][mm_plu]-new_r,2)
				     +0.5*pow(new_r-system[nn][mm_min],2)
				     //+0.5*pow(new_r,2);
						 +0.1*pow(new_r,4)
						 -1.0*pow(new_r,2);

				if(randomDouble(0,1)<exp(s_old-s_new)){ // metropolis acceptance
					system[nn][mm]=new_r;
	        accrate+=1.0/N_slices;
				}
			}
			if(i>500) h*=accrate/idrate;
			//cout << accrate << " " << h << endl;
			accrate=0.0;
		}

		// writing in the output file
		if(!(i%n_stride)){
			for(auto& particle : system){
				for(auto& pos : particle){
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
			cout << n/MCS << endl;
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
