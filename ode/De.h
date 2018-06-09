#pragma once
#include <cmath>
#include "Step.h"

class De
{
public:
	// constants

	// constructor and destructor
	De(void(*)(double, double[], double[]), unsigned char, double*, double, double, double, double, double*);
	~De();

	Step destep;

	// logical
	bool stiff;
	bool start;

	// variables
	unsigned char neqn;
	unsigned char nostep;
	unsigned char kle4;
	char isn;
	char isnold;
	char delsn;
	char iflag;

	short maxnum = 500;

	double u;
	double twou;
	double fouru;
	double relerr;
	double abserr;
	double eps;
	double releps;
	double abseps;
	double x;
	double xout;
	double t;
	double told;
	double tout;
	double tend;
	double h;
	double del;
	double absdel;

	// variable sized arrays;
	double* y;
	double* yp;
	double* yy;
	double* wt;
	double* ypout;
	double* yout;
	double* phi;

	// fixed size arrays
	double w[13];
	double g[13];
	double rho[13];

	// member functions
	void(*f)(double, double[], double[]);

	
	void step();

	void machine();

	void test_inputs();
	void setup();
	void first_step();
	void btout();
	void ctout();
	void work();
	void weights();
	void tolerances();
	void increment();
	void interp();

};

