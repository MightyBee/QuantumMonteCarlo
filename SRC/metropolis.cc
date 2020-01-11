#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <ctime>
#include "system.h"
using namespace std;


int main(int argc, char* argv[]){

	//############################## GENERATE RANDOM NUMBERS ##############################

	//Future code here. For now we just do
	srand(time(0));



	//############################## VILLARD LIBRARY ##############################

	string inputPath("../configuration.in"); // Fichier d'input par defaut
	if(argc>1) // Fichier d'input specifie par l'utilisateur ("./Exercice7 config_perso.in")
		inputPath = argv[1];

	ConfigFile configFile(inputPath); // Les parametres sont lus et stockes dans une "map" de strings.

	for(int i(2); i<argc; ++i) // Input complementaires ("./Exercice7 config_perso.in input_scan=[valeur]")
		configFile.process(argv[i]);


	//############################## READ PARAMETERS ##############################

	unsigned int N_sweeps(configFile.get<unsigned int>("N_sweeps"));		// number of Monte Carlo iterations (aka sweeps)
	double pos_min(configFile.get<double>("pos_min"));							// initial minimum position
	double pos_max(configFile.get<double>("pos_max"));							// initial maximal position
	double h(configFile.get<double>("h"));											// maximum uniform displacement of a point in the path
	double accrate(0.0);																	// ???
	double idrate(configFile.get<double>("idrate"));							// ???
	size_t n_stride(configFile.get<size_t>("n_stride"));						// output is written every n_stride iterations
	//Output file
	string output("../"+configFile.get<string>("output")+".out"); 					// output file
	ofstream fichier_output(output.c_str());
	fichier_output.precision(15);														// Precision

	System s(configFile);
	s.initialize(pos_min,pos_max);
	fichier_output << s << endl;


	//############################## METROPOLIS ALGORITHM ##############################
	//For every sweep...
	for(size_t i(0); i < N_sweeps; i++){
		//For every particle...
		// ??? should we directly make one 'for i=0:N_slices*N_part' ???
		for(size_t j(0); j < s.nb_part(); j++){
			for(size_t k(0); k < s.nb_slices(); k++){
				if(s.localMove(h)){
					accrate += 1.0/s.nb_slices();
				}
			}
			s.globalDisplacement(h);
			s.bissection(h);
			//if(i>500) h*=accrate/idrate;
			//cout << accrate << " " << h << endl;
			accrate = 0.0;
		}



		//############################## OUTPUT IN FILE ##############################
		if((i%n_stride) == 0){
			fichier_output << s << endl;
		}
	}
	fichier_output.close();



	//############################## STATS ABOUT VISITING ##############################
	/*//Impossible values so that we are sure to catch one point
	double visit_min(N_part * N_slices + 1.0);
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
	cout << "Most visited point, mean: " << visit_max << endl;*/

	//############################## END OF MAIN ##############################
	return 0;
}
