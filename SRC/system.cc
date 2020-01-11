#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include "system.h"
using namespace std;



System::System(const ConfigFile& configFile) :
	N_part(configFile.get<unsigned int>("N_part")),
	N_slices(configFile.get<unsigned int>("N_slices")),
	beta(configFile.get<double>("beta")),
	d_tau(beta/N_slices),
	mass(configFile.get<double>("mass")),
	omega(configFile.get<double>("frequency")),
	table(N_part, vector<double>(N_slices, 0.0)),
	mm(0), mm_plu(0), mm_min(0), nn(0),
	dis(0.0), s_old(0.0), s_new(0.0)
	{
		string type_V(configFile.get<string>("type_V"));
		if(type_V=="harmonic") ptr_V = new Potential_harm(configFile);
		else if(type_V=="double") ptr_V = new Potential_double(configFile);
		else if(type_V=="square") ptr_V = new Potential_square(configFile);
		else{
			cerr << "Please choose between ""harmonic"", ""double"" or ""square"" for ""type_V""." << endl;
		}
	}


void System::initialize(const double& pos_min, const double& pos_max){
	for(auto& particle : table){ // initialize random paths for each particles
		for(auto& pos : particle){
			pos = randomDouble(pos_min, pos_max);
		}
	}
}


ostream& System::write(ostream& output) const{
	for(const auto& particle : table){
		for(const auto& pos : particle){
			output << pos << " ";
		}
	}
	return output;
}


double System::kinetic(const int& particle, const int& bead, const int& bead_pm,
								const double& displacement) const{
	return 0.5*mass*pow(((table[particle][bead]+displacement)-table[particle][bead_pm])/d_tau,2);
}


bool System::localMove(const double& h){
	mm = rand()%N_slices; // random integer between 0 and N_slices-1
	mm_min = (mm + N_slices - 1)%N_slices; // mm-1 with periodic boundary condition
	mm_plu = (mm + 1)%N_slices; // mm+1 with periodic boundary condition
	nn = rand()%N_part; // random integer between 0 and N_part-1

	dis = h * randomDouble(-1, 1); // proposed new position

	// as we take the difference of new and old action S_new-S_old, we can
	// consider only the part of the action that is affected by the proposed new position
	s_old = kinetic(nn,mm,mm_plu) + kinetic(nn,mm,mm_min)
			+ (*ptr_V)(table[nn][mm]);
	s_new = kinetic(nn,mm,mm_plu,dis) + kinetic(nn,mm,mm_min,dis)
			+ (*ptr_V)(table[nn][mm]+dis);

	if(randomDouble(0,1) <= exp(-d_tau * (s_new - s_old))){ // metropolis acceptance
		table[nn][mm] += dis;		// update position with new one
		return true;
	}else{
		return false;
	}
}


bool System::globalDisplacement(const double& h){
	nn = rand()%N_part; // random integer between 0 and N_part-1
	dis = h * randomDouble(-1,1); // proposed displacement of the 'entire' particle nn

	// no relative move between the time slices --> only the potential action changes
	s_old=0.0;
	s_new=0.0;
	for(const auto& pos : table[nn]){
		s_old+=(*ptr_V)(pos);
		s_new+=(*ptr_V)(pos+dis);
	}

	if(randomDouble(0,1)<exp(-d_tau * (s_new - s_old))){ // metropolis acceptance
		for(auto& pos : table[nn]){
			pos+=dis;
		}
		return true;
	}else{
		return false;
	}
}



size_t System::nb_part() const {return N_part;}

size_t System::nb_slices() const {return N_slices;}


bool System::bissection(const double& h){
	mm = rand()%N_slices; // random integer between 0 and N_slices-1
	mm_min = (mm + N_slices - 1)%N_slices; // mm-1 with periodic boundary condition
	nn = rand()%N_part; // random integer between 0 and N_part-1
	dis = h * randomDouble(-1,1); // proposed displacement of the 'entire' particle nn
	size_t l(N_slices/3);

	// no relative move between the time slices --> only the potential action changes
	s_old=0.0;
	s_new=0.0;
	for(size_t i(0); i<l; i++){
		s_old+=(*ptr_V)(table[nn][(mm+i)%N_slices]);
		s_new+=(*ptr_V)(table[nn][(mm+i)%N_slices]+dis);
	}
	s_old += kinetic(nn,mm,mm_min)     + kinetic(nn,(mm+l-1)%N_slices,(mm+l)%N_slices);
	s_new += kinetic(nn,mm,mm_min,dis) + kinetic(nn,(mm+l-1)%N_slices,(mm+l)%N_slices,dis);

	if(randomDouble(0,1)<exp(-d_tau * (s_new - s_old))){ // metropolis acceptance
		for(size_t i(0); i<l; i++){
			table[nn][(mm+i)%N_slices]+=dis;
		}
		return true;
	}else{
		return false;
	}
}


ostream& operator<<(ostream& output, const System& s){
	return s.write(output);
}


//########################### FUNCTION DEFINITIONS #############################

double randomDouble(const double& min, const double& max){
	return (min + (max-min) * (double)rand()/RAND_MAX);
}
