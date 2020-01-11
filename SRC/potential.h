#include "instance.h"


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
