#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <ctime>
#include "ConfigFile.tcc"
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
double bissectionAction(const double& m, const double& omega, const double& d_tau, const size_t& start, const vector<double>& path, const double& displacement = 0.0);

// Potential function
double V(const double& x);		// V(x)
double dV(const double& x);	// V'(x)

//double QLagrangian(const vector<vector<double>>& pos, const double& d_tau);
//double diff_QLagrangian(const vector<vector<double>>& pos, const double& d_tau, const unsigned int& m, const unsigned int& n, const double& new_pos);


// Abstract class for potential
class Potential {
public:
	// pure virtual method => abstract class
	virtual double operator()(const double& x) const = 0; // return V at point x
};


// Class for a harmonic potential
class Potential_harm: public Potential {
public:
	Potential_harm(const ConfigFile& configFile);
	double operator()(const double& x) const;
private:
  double m, omega2;
};


// Class for a double well potential
class Potential_double: public Potential {
public:
	Potential_double(const ConfigFile& configFile);
	double operator()(const double& x) const;
private:
	double V0, x0;
};


// Class for a square potential (barrier for V0>0 and well for V0<0)
class Potential_square: public Potential {
public:
	Potential_square(const ConfigFile& configFile);
	double operator()(const double& x) const;
private:
	double V0, xc, L;
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
	bool localMove(const double& h);
	bool globalDisplacement(const double& h);
	bool bissection(const double& h);

private:
	// physical parameters
	unsigned int N_part;
	unsigned int N_slices;
	double beta;
	double d_tau;
	double mass;
	double omega;
	Potential* ptr_V;
	vector<vector<double>> table;
	// utilitary variables
	unsigned int mm, mm_plu, mm_min; // mm : time slice randomly selected during each iteration, mm_plu=mm+1, mm_min=mm-1;
	unsigned int nn; // particle randomly selected during each iteration
	double dis; // displacement proposed
	double s_old, s_new; // part of the action that is changed with the new position proposed (old and new values respectively)
};

ostream& operator<<(ostream& output, const System& s);


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

	unsigned int N_sweeps(configFile.get<unsigned int>("N_sweeps"));		// number of Monte Carlo iterations (aka sweeps)
	double pos_min(configFile.get<double>("pos_min"));							// initial minimum position
	double pos_max(configFile.get<double>("pos_max"));							// initial maximal position
	double h(configFile.get<double>("h"));											// maximum uniform displacement of a point in the path
	double accrate(0.0);																	// ???
	double idrate(configFile.get<double>("idrate"));							// ???
	size_t n_stride(configFile.get<size_t>("n_stride"));						// output is written every n_stride iterations
	//Output file
	string output(configFile.get<string>("output")+".out"); 					// output file
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

//########################### class 'Potential' methods DEFINITIONS #############################


//##### Potential_harm ######

Potential_harm::Potential_harm(const ConfigFile& configFile) :
	Potential(),
	m(configFile.get<double>("mass")),
	omega2(pow(configFile.get<double>("omega"),2))
	{}

double Potential_harm::operator()(const double& x) const {
	return m*omega2*x*x;
}


//##### Potential_double ######

Potential_double::Potential_double(const ConfigFile& configFile) :
	Potential(),
	V0(configFile.get<double>("V0")),
	x0(configFile.get<double>("x0"))
	{}

double Potential_double::operator()(const double& x) const {
	return V0*pow(pow(x/x0,2)-1,2);
}


//##### Potential_harm ######

Potential_square::Potential_square(const ConfigFile& configFile) :
	Potential(),
	V0(configFile.get<double>("V0")),
	xc(configFile.get<double>("xc")),
	L(configFile.get<double>("L"))
	{}

double Potential_square::operator()(const double& x) const {
	if(x<xc-0.5*L || x>xc+0.5*L){
		return 0;
	}else{
		return V0;
	}
}


//########################### class 'System' methods DEFINITIONS #############################

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
