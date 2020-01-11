#include <cmath>
#include "potential.h"
using namespace std;

//########################### class 'Potential' methods DEFINITIONS #############################



//##### Potential_harm ######

Potential_harm::Potential_harm(const ConfigFile& configFile) :
	Potential(),
	m(configFile.get<double>("mass")),
	omega2(pow(configFile.get<double>("omega"),2))
	{}

double Potential_harm::operator()(const double& x) const {
	return 0.5*m*omega2*x*x;
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
