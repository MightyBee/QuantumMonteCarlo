#define _USE_MATH_DEFINES
#include <vector>
#include <array>
#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <memory>
#include <ctime>
#include <cmath>
#include "ConfigFile.tcc"
using namespace std;

std::mt19937 rng(time(0));
//std::mt19937 rng(1);
///srand(time(0));

/* TODO :	finish swap move
				output with hamiltonian
*/

/*############################## NOTES ##############################
	- theoretically x_0, x_1, ... , x_(N_slices)					, (N_slices+1) points
	- we consider boundary conditions, x_0 = x_(N_slices)
	- hence we only consider x_0, x_1, ... , x_(N_slices-1)	, N_slices points
	- to get the "full picture" simply add one more point, equal to x_0
*/



//############################## HEADERS ##############################

// Generate a random (uniform) double between 'min' and 'max'
double randomDouble(const double& min=0.0, const double& max=1.0, const bool& closed=true);
// Generate a random double from a normal Cauchy distribution
double CauchyDistribution();
// Other distrubution ???

double GenerateDist(const double& h);

// 'true' modulo (in mathematic sense)
double mod(double a, double b){
	return a - b*floor(a/b);
}

// periodic condition for interval [-L/2,L/2]
double periodicBC(double x, double L){
	return mod(x+0.5*L,L)-0.5*L;
}



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
	double V0, x0, L;
};


// Class for a sinusoidal potential
class PotExt_sin: public Potential_ext {
public:
	PotExt_sin(const ConfigFile& configFile);
	double operator()(const double& x) const;
private:
	double V0, L;
};

class PotExt_LJ: public Potential_ext {
public:
	PotExt_LJ(const ConfigFile& configFile);
	double operator()(const double& x) const;
private:
	double V0, x0;
};

// Class for a sinusoidal potential
class PotExt_OHbonds: public Potential_ext {
public:
	PotExt_OHbonds(const ConfigFile& configFile);
	double Vmorse(const double& x) const;
	double operator()(const double& x) const;
private:
	double D, a, r0, delta1, b, R1;
	double R, DELTA;
};





 ////CAREFULLL HAVE TO TAKE BCs INTO ACCOUNT WHEN COMPUTAING INTERNAL POTENTIAL
// Abstract class for internal potential
class Potential_int {
public:
	// pure virtual method => abstract class
	virtual double operator()(const double& x1, const double& x2, const size_t& i, const size_t& j) const = 0; // return V at point x
};
 ////CAREFULLL HAVE TO TAKE BCs INTO ACCOUNT WHEN COMPUTAING INTERNAL POTENTIAL

// Class for a null internal potential (no interactions between particles)
class PotInt_null: public Potential_int {
public:
	double operator()(const double& x1, const double& x2, const size_t& i, const size_t& j) const {return 0.0;}
};


// Class for a harmonic potential between two particles
class PotInt_harm: public Potential_int {
public:
	PotInt_harm(const ConfigFile& configFile);
	double operator()(const double& x1, const double& x2, const size_t& i, const size_t& j) const;
private:
	double k, l0;
};


class PotInt_threeBonds: public Potential_int {
public:
	PotInt_threeBonds(const ConfigFile& configFile);
	double Vmorse(const double& r) const;
	double Vcoulomb(const double& r) const;
	double operator()(const double& x1, const double& x2, const size_t& i, const size_t& j) const;
private:
	double D, a, r0, A;
};



class System {
public:
	System(const ConfigFile& configFile);
	void initialize(const double& pos_min, const double& pos_max);
	void write_potExt(const string& output);
	size_t nb_part() const {return N_part;}
	size_t nb_slices() const {return N_slices;}
	vector<vector<int>> get_visits() const {return verif;}
	ostream& write(ostream& output) const;
	double kinetic(const int& particle, const int& bead, const int& bead_pm,
		 				const double& displacement=0.0) const;
	double energy();
	double get_H(){return H;}

	bool metropolisAcceptance();

	// different moves possible
	bool localMove(const double& h);
	bool globalDisplacement(const double& h);
	bool bissection(const double& h, const double& sRel);
	bool swap();
	bool inverse();
	bool symmetryCM();

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
	unique_ptr<Potential_ext> ptr_Vext;
	unique_ptr<Potential_int> ptr_Vint;
	vector<vector<double>> table;
	// utilitary variables
	unsigned int mm, mm_plu, mm_min; // mm : time slice randomly selected during each iteration, mm_plu=mm+1, mm_min=mm-1;
	unsigned int nn; // particle randomly selected during each iteration
	double dis; // displacement proposed
	double s_old, s_new; // part of the action that is changed with the new position proposed (old and new values respectively)
	double H;
	vector<double> energies_psi;
	vector<vector<int>> verif;
};

ostream& operator<<(ostream& output, const System& s);

void System::measure_energy(){
	double temp_energy(0);
	for(size_t i(0); i < table[0].size(); i++){
		temp_energy += pow(table[0][i], 2);
	}
	energies_psi.push_back(pow(omega, 2) * mass[0] * temp_energy/N_slices);
}

void System::average_energy(){
	ofstream fichier_output;
	fichier_output.open("energies.out");
	fichier_output.precision(15);

	double temp_energy(0), temp_error(0);

	cout << "Finally, with d_tau = " << d_tau << endl;

	for(size_t i(0); i < energies_psi.size(); i++){
		fichier_output << energies_psi[i] << endl;
		temp_energy += energies_psi[i];
		temp_error += pow(energies_psi[i], 2);
	}
	temp_energy = temp_energy/energies_psi.size();
	temp_error = sqrt((temp_error/energies_psi.size() - pow(temp_energy, 2))/energies_psi.size());

	cout << "PSI: " << temp_energy << " +- " << temp_error << endl;
	cout << "Theory: " << (10 * 0.5 * omega * 1.0545718) << ",e-20" << endl;

	fichier_output.close();
}


//################################################################################################//
//############################################  MAIN  ############################################//
//################################################################################################//


int main(int argc, char* argv[]){


	//############################## VILLARD LIBRARY ##############################

	string inputPath("configuration.in"); // Fichier d'input par defaut
	if(argc>1) // Fichier d'input specifie par l'utilisateur ("./Exercice7 config_perso.in")
		inputPath = argv[1];

	ConfigFile configFile(inputPath); // Les parametres sont lus et stockes dans une "map" de strings.

	for(int i(2); i<argc; ++i) // Input complementaires ("./Exercice7 config_perso.in input_scan=[valeur]")
		configFile.process(argv[i]);


	//############################## READ PARAMETERS ##############################

	unsigned int N_sweeps(configFile.get<unsigned int>("N_sweeps"));		// number of Monte Carlo iterations (aka sweeps)
	unsigned int N_thermalisation(configFile.get<unsigned int>("N_thermal"));	//How many steps should we wait before measuring stuff
	double pos_min(configFile.get<double>("pos_min"));							// initial minimum position
	double pos_max(configFile.get<double>("pos_max"));							// initial maximal position
	vector<double> h(3,configFile.get<double>("h"));							// maximum uniform displacement of a point in the path
	double p_loc(configFile.get<double>("p_local"));
	double p_dsp(configFile.get<double>("p_displacement"));
	double p_bis(configFile.get<double>("p_bissection"));
	double s_bis(configFile.get<double>("s_bissection"));
	double p_swap(configFile.get<double>("p_swap"));
	double p_inv(configFile.get<double>("p_inverse"));
	double p_sym(configFile.get<double>("p_symmetryCM"));
	vector<unsigned int> NbTries(6,0);												// numbers of tries for [0] local move [1] displacement and [2] bissection
	vector<double>			accrate(6,0.0);											// acceptance rates for the three moves
	double tmp_accrate(0.0);															// "instantaneous" acceptance rate for local moves
	double idrate(configFile.get<double>("idrate"));							// ideal acceptance rate
	size_t n_stride(configFile.get<size_t>("n_stride"));						// output is written every n_stride iterations
	//Output file
	string output(configFile.get<string>("output"));		 					// output file
	string output_pos(output+"_pos.out");
	ofstream fichier_output(output_pos.c_str());
	fichier_output.precision(15);														// Precision

	string output_energy(output+"_nrg.out");
	ofstream fichier_energy(output_energy.c_str());
	fichier_energy.precision(15);														// Precision

	System s(configFile);
	s.initialize(pos_min,pos_max);
	s.write_potExt(output);
	fichier_output << s << endl;
	fichier_energy << s.get_H() << " " << s.energy() << endl;



	//############################## METROPOLIS ALGORITHM ##############################
	double last_measured_time(time(0));

	//For every sweep...
	for(size_t i(0); i < N_sweeps; i++){

		if(time(0) - last_measured_time >= 5){
			last_measured_time = time(0);
			cout << floor((double)i/N_sweeps*100) << " %" << endl;
		}

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
			// swap
			if(NbTries[3]*1.0/(i*s.nb_part()+j+1) < p_swap){
				NbTries[3]++;
				if(s.swap()){
					accrate[3]++;
				}
			}
			// inverse
			if(NbTries[4]*1.0/(i*s.nb_part()+j+1) < p_inv){
				NbTries[4]++;
				if(s.inverse()){
					accrate[4]++;
				}
			}
			// symmetryCM
			if(NbTries[5]*1.0/(i*s.nb_part()+j+1) < p_sym){
				NbTries[5]++;
				if(s.symmetryCM()){
					accrate[5]++;
				}
			}
		}

		//cout << accrate[0] << ", " << h;
		//cout << " " << h << endl;

		//############################## OUTPUT IN FILE ##############################
		if((i%n_stride) == 0){
			fichier_output << s << endl;
			fichier_energy << s.get_H() << " " << s.energy() << endl; //" " << s.energy() << endl;

			//Energy measurement
			if(i >= N_thermalisation){
				s.measure_energy();
			}
		}
	}
	fichier_output.close();
	fichier_energy.close();

	//Statistics
	string output_stat(output+"_stat.out");
	fichier_output.open(output_stat.c_str());
	fichier_output.precision(15);
	for(size_t i(0); i<accrate.size(); i++){
		accrate[i]/=NbTries[i];
		fichier_output << NbTries[i] << " " << accrate[i] << endl;
	}

	for(const auto& part : s.get_visits()){
		for(const auto& v : part){
			fichier_output << v/(NbTries[0]*1.0/(s.nb_part()*s.nb_slices())) << " ";
		}
		fichier_output << endl;
	}

	fichier_output.close();

	//Energy
	s.average_energy();

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
	x0(configFile.get<double>("x0")),
	L(configFile.get<double>("L"))
	{}

double PotExt_square::operator()(const double& x) const {
	if(abs(x - x0) < L/2){
		return V0;
	}else{
		return 0;
	}
}


//##### PotExt_sin ######

PotExt_sin::PotExt_sin(const ConfigFile& configFile) :
	Potential_ext(),
	V0(configFile.get<double>("V0")),
	L(configFile.get<double>("L"))
	{}

double PotExt_sin::operator()(const double& x) const {
	return 0.5*V0*(1-cos(2*M_PI/L*x));
}


PotExt_LJ::PotExt_LJ(const ConfigFile& configFile) :
	Potential_ext(),
	V0(configFile.get<double>("V0")),
	x0(configFile.get<double>("x0"))
	{}

double PotExt_LJ::operator()(const double& x) const {
	return 4 * V0 * (pow(x/x0, 12) - pow(x/x0, 6) );
}

PotExt_OHbonds::PotExt_OHbonds(const ConfigFile& configFile) :
	Potential_ext(),
	D(83.402), a(2.2), r0(0.96), delta1(0.4*D), b(2.2), R1(2*r0+1/a),
	R(configFile.get<double>("R")),
	DELTA(delta1*exp(-b*(R-R1)))
	{}
double PotExt_OHbonds::Vmorse(const double& x) const{
	return D*(exp(-2*a*(x-r0))-2*exp(-a*(x-r0)));
}
double PotExt_OHbonds::operator()(const double& x) const{
	if(abs(x) < 8){
		return 0.5*(Vmorse(R/2+x)+Vmorse(R/2-x) - sqrt(pow(Vmorse(R/2+x)-Vmorse(R/2-x),2)+4*DELTA*DELTA));
	}else{
		return 0.0;
	}
}

//##### Potential_rel ######

PotInt_harm::PotInt_harm(const ConfigFile& configFile) :
	Potential_int(),
	k(configFile.get<double>("k")),
	l0(configFile.get<double>("l0"))
	{}

 ////CAREFULLL HAVE TO TAKE BCs INTO ACCOUNT WHEN COMPUTAING INTERNAL POTENTIAL
double PotInt_harm::operator()(const double& x1, const double& x2, const size_t& i, const size_t& j) const {
	return 0.5*k*pow(abs(x2-x1)-l0,2);
}


PotInt_threeBonds::PotInt_threeBonds(const ConfigFile& configFile) :
	Potential_int(),
	D(83.402), a(2.2), r0(0.96),
	A(pow(configFile.get<double>("Q"),2)*230.708)
	{}

double PotInt_threeBonds::Vmorse(const double& r) const{
	return D*(exp(-2*a*(r-r0))-2*exp(-a*(r-r0)));
}

double PotInt_threeBonds::Vcoulomb(const double& r) const {
	return A/r;
}

 ////CAREFULLL HAVE TO TAKE BCs INTO ACCOUNT WHEN COMPUTAING INTERNAL POTENTIAL
double PotInt_threeBonds::operator()(const double& x1, const double& x2, const size_t& i, const size_t& j) const {
	if(i==2 or j==2){
		return Vmorse(abs(x1-x2));
	}
	else{
		return Vcoulomb(abs(x1-x2));
	}
}


//########################### class 'System' methods DEFINITIONS #############################

System::System(const ConfigFile& configFile) :
	N_part(configFile.get<unsigned int>("N_part")),
	N_slices(configFile.get<unsigned int>("N_slices")),
	beta(configFile.get<double>("beta")),
	d_tau(beta/N_slices),
	mass(N_part,0.0),// configFile.get<double>("mass")),
	omega(configFile.get<double>("omega")),
	table(N_part, vector<double>(N_slices, 0.0)),
	old_table(table),
	mm(0), mm_plu(0), mm_min(0), nn(0),
	dis(0.0), s_old(0.0), s_new(0.0), H(0.0),
	verif(N_part, vector<int>(N_slices, 0))
	{
		for(unsigned int i(0); i<N_part; i++){
			mass[i]=configFile.get<double>("m"+to_string(i+1));
		}

		string V_ext(configFile.get<string>("V_ext"));
		if(V_ext=="null") ptr_Vext = move(unique_ptr<Potential_ext>(new PotExt_null()));
		else if(V_ext=="harmonic") ptr_Vext = move(unique_ptr<Potential_ext>(new PotExt_harm(configFile)));
		else if(V_ext=="double") ptr_Vext = move(unique_ptr<Potential_ext>(new PotExt_double(configFile)));
		else if(V_ext=="square") ptr_Vext = move(unique_ptr<Potential_ext>(new PotExt_square(configFile)));
		else if(V_ext=="sin") ptr_Vext = move(unique_ptr<Potential_ext>(new PotExt_sin(configFile)));
		else if(V_ext=="LJ") ptr_Vext = move(unique_ptr<Potential_ext>(new PotExt_LJ(configFile)));
		else if(V_ext=="OHbonds") ptr_Vext = move(unique_ptr<Potential_ext>(new PotExt_OHbonds(configFile)));
		else{
			cerr << "Please choose a valid potential." << endl;
		}
		string V_int(configFile.get<string>("V_int"));
		if(V_int=="null") ptr_Vint = move(unique_ptr<Potential_int>(new PotInt_null()));
		else if(V_int=="harmonic") ptr_Vint = move(unique_ptr<Potential_int>(new PotInt_harm(configFile)));
		else if(V_int=="threeBonds") ptr_Vint = move(unique_ptr<Potential_int>(new PotInt_threeBonds(configFile)));
		else{
			cerr << "Please choose a valid potential." << endl;
		}
	}


void System::initialize(const double& pos_min, const double& pos_max){
	for(auto& particle : table){ // initialize random paths for each particles
		for(auto& pos : particle){
			pos = randomDouble(pos_min, pos_max);
		}
	}
	H=energy();
}


void System::write_potExt(const string& output){
	string output_pot(output+"_pot.out");
	ofstream f_pot(output_pot.c_str());
	f_pot.precision(15);
	size_t N(10000);
	double x(0.0);
	for(size_t i(0); i<N; i++){
		double xi(-200.0), xf(200.0);
		x=xi+i*(xf-xi)/(N-1);
		f_pot << x << " " << (*ptr_Vext)(x) << endl;
	}
	f_pot.close();
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

double System::energy(){
	double E(0.0);
	for(size_t part(0); part<N_part; part++){
		for(size_t bead(0); bead<N_slices; bead++){
			E+=kinetic(part,bead,(bead+1)%N_slices);
			E+=(*ptr_Vext)(table[part][bead]);
			for(size_t part2(part+1); part2<N_part; part2++){
				E+=(*ptr_Vint)(table[part][bead],table[part2][bead],part,part2);
			}
		}
	}
	return E;
}



bool System::metropolisAcceptance(){
	//return (randomDouble(0,1) <= exp(-(d_tau * pow(10, -20-q))* (s_new - s_old)));
	//cout << (-(d_tau) * (s_new - s_old)) << endl;
	return ( randomDouble(0,1) <= exp(-(0.1*d_tau/1.0545718) * (s_new - s_old)) );
	//return (randomDouble(0,1) <= exp(-(s_new - s_old)));
}


bool System::localMove(const double& h){
	mm = rng()%N_slices; // random integer between 0 and N_slices-1
	mm_min = (mm + N_slices - 1)%N_slices; // mm-1 with periodic boundary condition
	mm_plu = (mm + 1)%N_slices; // mm+1 with periodic boundary condition
	nn = rng()%N_part; // random integer between 0 and N_part-1

	verif[nn][mm]++;

	dis=GenerateDist(h);

	// as we take the difference of new and old action S_new-S_old, we can
	// consider only the part of the action that is affected by the proposed new position
	s_old = kinetic(nn,mm,mm_plu) + kinetic(nn,mm,mm_min)
			+ (*ptr_Vext)(table[nn][mm]);
	s_new = kinetic(nn,mm,mm_plu,dis) + kinetic(nn,mm,mm_min,dis)
			+ (*ptr_Vext)(table[nn][mm]+dis);

	if(N_part>1){
		for(size_t i(0); i<N_part; i++){
			if(i!=nn){
				s_old+=(*ptr_Vint)(table[i][mm],table[nn][mm],i,nn);
				s_new+=(*ptr_Vint)(table[i][mm],table[nn][mm]+dis,i,nn);
			}
		}
	}


	if(metropolisAcceptance()){ // metropolis acceptance
		table[nn][mm] += dis;		// update position with new one
		H+=s_new-s_old;
		return true;
	}else{
		return false;
	}
}


bool System::globalDisplacement(const double& h){
	nn = rng()%N_part; // random integer between 0 and N_part-1
	dis=GenerateDist(h);

	// no relative move between the time slices --> only the potential action changes
	s_old=0.0;
	s_new=0.0;
	for(size_t j(0); j<N_slices; j++){
		s_old+=(*ptr_Vext)(table[nn][j]);
		s_new+=(*ptr_Vext)(table[nn][j]+dis);
		if(N_part>1){
			for(size_t i(0); i<N_part; i++){
				if(i!=nn){
					s_old+=(*ptr_Vint)(table[i][j],table[nn][j],i,nn);
					s_new+=(*ptr_Vint)(table[i][j],table[nn][j]+dis,i,nn);
				}
			}
		}
	}

	if(metropolisAcceptance()){ // metropolis acceptance
		for(auto& pos : table[nn]){
			pos+=dis;
		}
		H+=s_new-s_old;
		return true;
	}else{
		return false;
	}
}



bool System::bissection(const double& h, const double& sRel){
	mm = rng()%N_slices; // random integer between 0 and N_slices-1
	mm_min = (mm + N_slices - 1)%N_slices; // mm-1 with periodic boundary condition
	nn = rng()%N_part; // random integer between 0 and N_part-1
	dis=GenerateDist(h);
	size_t l(N_slices*sRel);


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
					s_old+=(*ptr_Vint)(table[i][ind_j],table[nn][ind_j],i,nn);
					s_new+=(*ptr_Vint)(table[i][ind_j],table[nn][ind_j]+dis,i,nn);
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
		H+=s_new-s_old;
		return true;
	}else{
		return false;
	}
}



bool System::swap(){
	if (N_part>1){
		mm_min = rng()%N_part; // random integer between 0 and N_part-1
		mm_plu = (mm_min+1+rng()%(N_part-1))%N_part; // another, but different, random integer between 0 and N_part-1
		mm = rng()%N_slices; // random integer between 0 and N_slices-1 (bead where the swap starts)
		nn = rng()%(N_slices-1)+1; //length of the swap (nb of slices swapped)

		// no relative move between the time slices --> only the potential action changes
		s_old=0.0;
		s_new=0.0;
		int ind_j(0), ind_j_pm(0);

		ind_j=mm;
		ind_j_pm=(mm+N_slices-1)%N_slices;
		s_old += kinetic(mm_min,ind_j,ind_j_pm) + kinetic(mm_plu,ind_j,ind_j_pm);
		s_new += kinetic(mm_min,ind_j,ind_j_pm,table[mm_plu][mm]-table[mm_min][mm]) + kinetic(mm_plu,ind_j,ind_j_pm,table[mm_min][mm]-table[mm_plu][mm]);

		for(size_t j(0); j<nn; j++){
			ind_j=(mm+j)%N_slices;
			ind_j_pm=(ind_j+N_slices-1)%N_slices;
			if(mass[mm_min]!=mass[mm_plu]){
				for(size_t i(0); i<N_part; i++){
					if(i!=mm_min and i!=mm_plu){
						// change for particle(mm_min)
						s_old+=(*ptr_Vint)(table[i][ind_j],table[mm_min][ind_j],i,mm_min);
						s_new+=(*ptr_Vint)(table[i][ind_j],table[mm_plu][ind_j],i,mm_min);
						// change for particle(mm_plu)
						s_old+=(*ptr_Vint)(table[i][ind_j],table[mm_plu][ind_j],i,mm_plu);
						s_new+=(*ptr_Vint)(table[i][ind_j],table[mm_min][ind_j],i,mm_plu);
					}
				}
			}
			if(j){ // in swapped part of paths : K1_new = m1*(K2_old/m2), K2_new=m2/m1*K1_old
				s_old += kinetic(mm_min,ind_j,ind_j_pm) + kinetic(mm_plu,ind_j,ind_j_pm);
				s_new += mass[mm_min]/mass[mm_plu]*kinetic(mm_plu,ind_j,ind_j_pm) + mass[mm_plu]/mass[mm_min]*kinetic(mm_min,ind_j,ind_j_pm);
			}
		}
		ind_j=(mm+nn-1)%N_slices;
		ind_j_pm=(mm+nn)%N_slices;
		s_old += kinetic(mm_min,ind_j,ind_j_pm) + kinetic(mm_plu,ind_j,ind_j_pm);
		s_new += kinetic(mm_min,ind_j,ind_j_pm,table[mm_plu][ind_j]-table[mm_min][ind_j]) + kinetic(mm_plu,ind_j,ind_j_pm,table[mm_min][ind_j]-table[mm_plu][ind_j]);


		if(metropolisAcceptance()){ // metropolis acceptance
			double tmp(0.0);
			for(size_t j(0); j<nn; j++){
				ind_j=(mm+j)%N_slices;
				tmp=table[mm_min][ind_j];
				table[mm_min][ind_j]=table[mm_plu][ind_j];
				table[mm_plu][ind_j]=tmp;
			}

			H+=s_new-s_old;
			return true;
		}
	}
	return false;
}

bool System::inverse(){
	nn = rng()%N_part; // random integer between 0 and N_part-1

	// no relative move between the time slices --> only the potential action changes
	s_old=0.0;
	s_new=0.0;
	for(size_t j(0); j<N_slices; j++){
		s_old+=(*ptr_Vext)( table[nn][j]);
		s_new+=(*ptr_Vext)(-table[nn][j]);
		if(N_part>1){
			for(size_t i(0); i<N_part; i++){
				if(i!=nn){
					s_old+=(*ptr_Vint)(table[i][j], table[nn][j],i,nn);
					s_new+=(*ptr_Vint)(table[i][j],-table[nn][j],i,nn);
				}
			}
		}
	}

	if(metropolisAcceptance()){ // metropolis acceptance
		for(auto& pos : table[nn]){
			pos*=-1;
		}
		H+=s_new-s_old;
		return true;
	}else{
		return false;
	}
}


bool System::symmetryCM(){
	nn = rng()%N_part; // random integer between 0 and N_part-1
	dis= 0;
	for(const auto& pos : table[nn]){
		dis+=pos;
	}
	dis*=-2.0/table[nn].size();

	// no relative move between the time slices --> only the potential action changes
	s_old=0.0;
	s_new=0.0;
	for(size_t j(0); j<N_slices; j++){
		s_old+=(*ptr_Vext)(table[nn][j]);
		s_new+=(*ptr_Vext)(table[nn][j]+dis);
		if(N_part>1){
			for(size_t i(0); i<N_part; i++){
				if(i!=nn){
					s_old+=(*ptr_Vint)(table[i][j],table[nn][j],i,nn);
					s_new+=(*ptr_Vint)(table[i][j],table[nn][j]+dis,i,nn);
				}
			}
		}
	}

	if(metropolisAcceptance()){ // metropolis acceptance
		for(auto& pos : table[nn]){
			pos+=dis;
		}
		H+=s_new-s_old;
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
	if(closed) return (min + (max-min) * (double)rng()/rng.max());
	else return (min + (max-min) * ((double)rng()+0.5)/(rng.max()+1.0));
}

double CauchyDistribution(){
	return tan(M_PI*(randomDouble(-0.5,0.5,false)));
}

double GenerateDist(const double& h){
	if(rng()%2){
		return h * randomDouble(-1.0,1.0); // proposed displacement
	}else{
		return h * CauchyDistribution(); // proposed displacement
	}
}
