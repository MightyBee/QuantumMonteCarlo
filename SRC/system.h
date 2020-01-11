#include <vector>
#include <iostream>
#include "potential.h"

class System {
public:
	System(const ConfigFile& configFile);
	void initialize(const double& pos_min, const double& pos_max);
	size_t nb_part() const;
	size_t nb_slices() const;
	std::ostream& write(std::ostream& output) const;
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
	std::vector<std::vector<double>> table;
	// utilitary variables
	unsigned int mm, mm_plu, mm_min; // mm : time slice randomly selected during each iteration, mm_plu=mm+1, mm_min=mm-1;
	unsigned int nn; // particle randomly selected during each iteration
	double dis; // displacement proposed
	double s_old, s_new; // part of the action that is changed with the new position proposed (old and new values respectively)
};

std::ostream& operator<<(std::ostream& output, const System& s);

// Generate a random (uniform) double between 'min' and 'max'
double randomDouble(const double& min, const double& max);

//
