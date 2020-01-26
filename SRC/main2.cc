#define _USE_MATH_DEFINES
#include <vector>
#include <array>
#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <ctime>
#include <cmath>
#include "ConfigFile.tcc"
using namespace std;



/*############################## NOTES ##############################
	- theoretically x_0, x_1, ... , x_(N_slices)					, (N_slices+1) points
	- we consider boundary conditions, x_0 = x_(N_slices)
	- hence we only consider x_0, x_1, ... , x_(N_slices-1)	, N_slices points
	- to get the "full picture" simply add one more point, equal to x_0
*/



//############################## HEADERS ##############################

// Generate a random (uniform) double between 'min' and 'max'
double randomDouble(const double& min, const double& max, const bool& closed=true);
// Generate a random double from a normal Cauchy distribution
double CauchyDistribution();


//double QLagrangian(const vector<vector<double>>& pos, const double& d_tau);
//double diff_QLagrangian(const vector<vector<double>>& pos, const double& d_tau, const unsigned int& m, const unsigned int& n, const double& new_pos);


// Abstract class for external potential
class Potential_ext {
public:
	// pure virtual method => abstract class
	virtual double operator()(const double& x) const = 0; // return V at point x
};


// Class for a null potential
class PotExt_null: public Potential_ext {
public:
	double operator()(const double& x) const {return 0.0;}
};


// Class for a harmonic potential
class PotExt_harm: public Potential_ext {
public:
	PotExt_harm(const ConfigFile& configFile);
	double operator()(const double& x) const;
private:
  double m, omega2;
};


// Class for a double well potential
class PotExt_double: public Potential_ext {
public:
	PotExt_double(const ConfigFile& configFile);
	double operator()(const double& x) const;
private:
	double V0, x0;
};


// Class for a square potential (barrier for V0>0 and well for V0<0)
class PotExt_square: public Potential_ext {
public:
	PotExt_square(const ConfigFile& configFile);
	double operator()(const double& x) const;
private:
	double V0, xc, L;
};


// Abstract class for internal potential
class Potential_int {
public:
	// pure virtual method => abstract class
	virtual double operator()(const double& x1, const double& x2) const = 0; // return V at point x
};


// Class for a null internal potential (no interactions between particles)
class PotInt_null: public Potential_int {
public:
	double operator()(const double& x1, const double& x2) const {return 0.0;}
};


// Class for a harmonic potential between two particles
class PotInt_harm: public Potential_int {
public:
	PotInt_harm(const ConfigFile& configFile);
	double operator()(const double& x1, const double& x2) const;
private:
	double k, l0;
};



class System {
public:
	System(const ConfigFile& configFile);
	void initialize(const double& pos_min, const double& pos_max);
	size_t nb_part() const {return N_part;}
	size_t nb_slices() const {return N_slices;}
	ostream& write(ostream& output) const;
	double kinetic(const int& particle, const int& bead, const int& bead_pm,
		 				const double& displacement=0.0) const;
	bool metropolisAcceptance();
	bool localMove(const double& h);
	bool globalDisplacement(const double& h);
	bool bissection(const double& h, const double& sRel);

	void measure_energy();
	void average_energy();

private:
	// physical parameters
	unsigned int N_part;
	unsigned int N_slices;
	double beta;
	double d_tau;
	int q;
	vector<double> mass;
	double omega;
	Potential_ext* ptr_Vext;
	Potential_int* ptr_Vint;
	vector<vector<double>> table;
	// utilitary variables
	unsigned int mm, mm_plu, mm_min; // mm : time slice randomly selected during each iteration, mm_plu=mm+1, mm_min=mm-1;
	unsigned int nn; // particle randomly selected during each iteration
	double dis; // displacement proposed
	double s_old, s_new; // part of the action that is changed with the new position proposed (old and new values respectively)
	vector<double> energies_psi;
	vector<double> energies_theo;
};

ostream& operator<<(ostream& output, const System& s);

void System::measure_energy(){
	double temp_energy1(0), temp_energy2(0);
	for(size_t i(0); i < table[0].size(); i++){
		temp_energy1 += pow(table[0][i], 2);
	}
	for(size_t i(0); i < table[0].size(); i++){
		temp_energy2 += (pow(table[0][i-1], 2) + pow(table[0][i], 2));
	}
	energies_psi.push_back(pow(omega, 2) * mass[0] * temp_energy1/N_slices);
	energies_theo.push_back(0.5 * pow(omega, 2) * mass[0] * temp_energy2/N_slices);
}

void System::average_energy(){
	double temp_energy(0);

	cout << "Finally, with d_tau = " << d_tau << " and m = " << mass[0] << endl;

	for(size_t i(0); i < energies_psi.size(); i++){		// !!!!! at the beginning i=1 and not i=0
		temp_energy += energies_psi[i];
	}

	cout << "PSI: " << (temp_energy/energies_psi.size() * pow(10, -20)) << endl;

	temp_energy = 0;
	for(size_t i(0); i < energies_theo.size(); i++){		// !!!!! at the beginning i=1 and not i=0
		temp_energy += energies_theo[i];
	}

	cout << "Theo: " << (temp_energy/energies_theo.size() * pow(10, -20)) << endl;
	cout << "Theory: " << (omega/d_tau * 0.5 * 1.0545718 * pow(10, q-34)) << endl;
}


//################################################################################################//
//############################################  MAIN  ############################################//
//################################################################################################//


int main(int argc, char* argv[]){

	//############################## GENERATE RANDOM NUMBERS ##############################

	//Future code here. For now we just do
	//srand(time(0));
	srand(100);



	//############################## VILLARD LIBRARY ##############################

	string inputPath("configuration.in"); // Fichier d'input par defaut
	if(argc>1) // Fichier d'input specifie par l'utilisateur ("./Exercice7 config_perso.in")
		inputPath = argv[1];

	ConfigFile configFile(inputPath); // Les parametres sont lus et stockes dans une "map" de strings.

	for(int i(2); i<argc; ++i) // Input complementaires ("./Exercice7 config_perso.in input_scan=[valeur]")
		configFile.process(argv[i]);


	//############################## READ PARAMETERS ##############################

	unsigned int N_sweeps(configFile.get<unsigned int>("N_sweeps"));		// number of Monte Carlo iterations (aka sweeps)
	unsigned int N_thermalisation(configFile.get<unsigned int>("N_thermalisation"));	//How many steps should we wait before measuring stuff
	double pos_min(configFile.get<double>("pos_min"));							// initial minimum position
	double pos_max(configFile.get<double>("pos_max"));							// initial maximal position
	vector<double> h(3,configFile.get<double>("h"));							// maximum uniform displacement of a point in the path
	double p_loc(configFile.get<double>("p_local"));
	double p_dsp(configFile.get<double>("p_displacement"));
	double p_bis(configFile.get<double>("p_bissection"));
	double s_bis(configFile.get<double>("s_bissection"));
	vector<unsigned int> NbTries(3,0);												// numbers of tries for [0] local move [1] displacement and [2] bissection
	vector<double>			accrate({0.0,0.0,0.0});									// acceptance rates for the three moves
	double tmp_accrate(0.0);															// "instantaneous" acceptance rate for local moves
	double idrate(configFile.get<double>("idrate"));							// ideal acceptance rate
	size_t n_stride(configFile.get<size_t>("n_stride"));						// output is written every n_stride iterations
	//Output file
	string output(configFile.get<string>("output"));		 					// output file
	string output_pos(output+"_pos.out");
	ofstream fichier_output(output_pos.c_str());
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
			// local move
			if(NbTries[0]*1.0/((i*s.nb_part()+j+1)*s.nb_slices()) < p_loc){
				tmp_accrate=0.0;
				for(size_t k(0); k < s.nb_slices(); k++){
					NbTries[0]++;
					if(s.localMove(h[0])){
						accrate[0]++;
						tmp_accrate++;
					}
				}
				tmp_accrate/=s.nb_slices();
				h[0]*=tmp_accrate/idrate;
			}
			// global displacement
			if(NbTries[1]*1.0/(i*s.nb_part()+j+1) < p_dsp){
				NbTries[1]++;
				if(s.globalDisplacement(h[1])){
					accrate[1]++;
				}
			}
			// bissection
			if(NbTries[2]*1.0/(i*s.nb_part()+j+1) < p_bis){
				NbTries[2]++;
				if(s.bissection(h[2], s_bis)){
					accrate[2]++;
				}
			}
		}
		//############################## OUTPUT IN FILE ##############################
		if((i%n_stride) == 0){
			fichier_output << s << endl;

			//Energy measurement
			if(i >= N_thermalisation){
				s.measure_energy();
			}
		}
	}
	fichier_output.close();

	//Statistics
	string output_stat(output+"_stat.out");
	fichier_output.open(output_stat.c_str());
	fichier_output.precision(15);
	for(size_t i(0); i<accrate.size(); i++){
		accrate[i]/=NbTries[i];
		fichier_output << NbTries[i] << " " << accrate[i] << endl;
	}
	fichier_output.close();

	//Energy
	s.average_energy();



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


//################################################################################################//
//#########################################  END OF MAIN  ########################################//
//################################################################################################//



//########################### class 'Potential_ext' methods DEFINITIONS #############################


//##### PotExt_harm ######

PotExt_harm::PotExt_harm(const ConfigFile& configFile) :
	Potential_ext(),
	m(configFile.get<double>("mass")),
	omega2(pow(configFile.get<double>("omega"),2))
	{}

double PotExt_harm::operator()(const double& x) const {
	return 0.5 * m * pow(x, 2) * omega2;
}


//##### PotExt_double ######

PotExt_double::PotExt_double(const ConfigFile& configFile) :
	Potential_ext(),
	V0(configFile.get<double>("V0")),
	x0(configFile.get<double>("x0"))
	{}

double PotExt_double::operator()(const double& x) const {
	return V0*pow(pow(x/x0,2)-1,2);
}


//##### PotExt_square ######

PotExt_square::PotExt_square(const ConfigFile& configFile) :
	Potential_ext(),
	V0(configFile.get<double>("V0")),
	xc(configFile.get<double>("xc")),
	L(configFile.get<double>("L"))
	{}

double PotExt_square::operator()(const double& x) const {
	if(x<xc-0.5*L || x>xc+0.5*L){
		return 0;
	}else{
		return V0;
	}
}


//##### Potential_rel ######

PotInt_harm::PotInt_harm(const ConfigFile& configFile) :
	Potential_int(),
	k(configFile.get<double>("k")),
	l0(configFile.get<double>("l0"))
	{}

double PotInt_harm::operator()(const double& x1, const double& x2) const {
	return 0.5*k*pow(abs(x2-x1)-l0,2);
}


//########################### class 'System' methods DEFINITIONS #############################

System::System(const ConfigFile& configFile) :
	N_part(configFile.get<unsigned int>("N_part")),
	N_slices(configFile.get<unsigned int>("N_slices")),
	beta(configFile.get<double>("beta")),
	d_tau(beta/N_slices),
	q(configFile.get<int>("q")),
	mass(N_part, configFile.get<double>("mass")),
	omega(configFile.get<double>("omega")),
	table(N_part, vector<double>(N_slices, 0.0)),
	mm(0), mm_plu(0), mm_min(0), nn(0),
	dis(0.0), s_old(0.0), s_new(0.0)
	{
		for(int i(0); i<N_part; i++){
			mass[i]=configFile.get<double>("m"+to_string(i+1));
		}
		string V_ext(configFile.get<string>("V_ext"));
		if(V_ext=="null") ptr_Vext = new PotExt_null();
		else if(V_ext=="harmonic") ptr_Vext = new PotExt_harm(configFile);
		else if(V_ext=="double") ptr_Vext = new PotExt_double(configFile);
		else if(V_ext=="square") ptr_Vext = new PotExt_square(configFile);
		else{
			cerr << "Please choose between ""null"", ""harmonic"", ""double"" or ""square"" for ""V_ext""." << endl;
		}
		string V_int(configFile.get<string>("V_int"));
		if(V_int=="null") ptr_Vint = new PotInt_null();
		else if(V_int=="harmonic") ptr_Vint = new PotInt_harm(configFile);
		else{
			cerr << "Please choose between ""null"", ""harmonic"" for ""V_int""." << endl;
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


double System::kinetic(const int& particle, const int& bead, const int& bead_pm, const double& displacement) const{
	return 0.5*mass[particle]*pow(((table[particle][bead]+displacement)-table[particle][bead_pm])/d_tau,2);
	//return 0.5*mass[particle]*pow((table[particle][bead]+displacement)-table[particle][bead_pm],2);
}



bool System::metropolisAcceptance(){
	//return (randomDouble(0,1) <= exp(-(d_tau * pow(10, -20-q))* (s_new - s_old)));
	//cout << (-(d_tau) * (s_new - s_old)) << endl;
	return (randomDouble(0,1) <= exp(-(d_tau) * (s_new - s_old)));
	//return (randomDouble(0,1) <= exp(-(s_new - s_old)));
}


bool System::localMove(const double& h){
	mm = rand()%N_slices; // random integer between 0 and N_slices-1
	mm_min = (mm + N_slices - 1)%N_slices; // mm-1 with periodic boundary condition
	mm_plu = (mm + 1)%N_slices; // mm+1 with periodic boundary condition
	nn = rand()%N_part; // random integer between 0 and N_part-1

	//dis = h * randomDouble(-1, 1); // proposed new position
	dis = h * CauchyDistribution(); // proposed new position (Cauchy distribution)

	// as we take the difference of new and old action S_new-S_old, we can
	// consider only the part of the action that is affected by the proposed new position
	s_old = kinetic(nn,mm,mm_plu) + kinetic(nn,mm,mm_min)
			+ (*ptr_Vext)(table[nn][mm]);
	s_new = kinetic(nn,mm,mm_plu,dis) + kinetic(nn,mm,mm_min,dis)
			+ (*ptr_Vext)(table[nn][mm]+dis);
	if(N_part>1){
		for(size_t i(0); i<N_part; i++){
			if(i!=nn){
				s_old+=(*ptr_Vint)(table[i][mm],table[nn][mm]);
				s_new+=(*ptr_Vint)(table[i][mm],table[nn][mm]+dis);
			}
		}
	}

	if(metropolisAcceptance()){ // metropolis acceptance
		table[nn][mm] += dis;		// update position with new one
		return true;
	}else{
		return false;
	}
}


bool System::globalDisplacement(const double& h){
	nn = rand()%N_part; // random integer between 0 and N_part-1
	//dis = h * randomDouble(-1,1); // proposed displacement of the 'entire' particle nn
	dis = h * CauchyDistribution();  // proposed displacement of the 'entire' particle nn

	// no relative move between the time slices --> only the potential action changes
	s_old=0.0;
	s_new=0.0;
	for(size_t j(0); j<N_slices; j++){
		s_old+=(*ptr_Vext)(table[nn][j]);
		s_new+=(*ptr_Vext)(table[nn][j]+dis);
		if(N_part>1){
			for(size_t i(0); i<N_part; i++){
				if(i!=nn){
					s_old+=(*ptr_Vint)(table[i][j],table[nn][j]);
					s_new+=(*ptr_Vint)(table[i][j],table[nn][j]+dis);
				}
			}
		}
	}

	if(metropolisAcceptance()){ // metropolis acceptance
		for(auto& pos : table[nn]){
			pos+=dis;
		}
		return true;
	}else{
		return false;
	}
}



bool System::bissection(const double& h, const double& sRel){
	mm = rand()%N_slices; // random integer between 0 and N_slices-1
	mm_min = (mm + N_slices - 1)%N_slices; // mm-1 with periodic boundary condition
	nn = rand()%N_part; // random integer between 0 and N_part-1
	//dis = h * randomDouble(-1,1); // proposed displacement
	dis = h * CauchyDistribution(); // proposed displacement
	size_t l(N_slices*sRel);

	// no relative move between the time slices --> only the potential action changes
	s_old=0.0;
	s_new=0.0;
	int ind_j(0);
	for(size_t j(0); j<l; j++){
		ind_j=(mm+j)%N_slices;
		s_old+=(*ptr_Vext)(table[nn][ind_j]);
		s_new+=(*ptr_Vext)(table[nn][ind_j]+dis);
		if(N_part>1){
			for(size_t i(0); i<N_part; i++){
				if(i!=nn){
					s_old+=(*ptr_Vint)(table[i][ind_j],table[nn][ind_j]);
					s_new+=(*ptr_Vint)(table[i][ind_j],table[nn][ind_j]+dis);
				}
			}
		}
	}
	s_old += kinetic(nn,mm,mm_min)     + kinetic(nn,(mm+l-1)%N_slices,(mm+l)%N_slices);
	s_new += kinetic(nn,mm,mm_min,dis) + kinetic(nn,(mm+l-1)%N_slices,(mm+l)%N_slices,dis);

	if(metropolisAcceptance()){ // metropolis acceptance
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

double randomDouble(const double& min, const double& max, const bool& closed){
	if(closed) return (min + (max-min) * (double)rand()/RAND_MAX);
	else return (min + (max-min) * ((double)rand()+0.5)/(RAND_MAX+1.0));
}

double CauchyDistribution(){
	return tan(M_PI*(randomDouble(-0.5,0.5,false)));
}
