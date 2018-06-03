#pragma once
#include <cmath>
#include "Step.h"

class De
{
public:
	// constants

	// constructor and destructor
	De();
	~De();

	Step destep;

	// logical
	bool stiff;
	bool start;

	// variables
	char neqn;
	char nostep;
	char kle4;
	char isn;
	char delsn;
	char iflag;

	double relerr;
	double abserr;
	double eps;
	double releps;
	double abseps;
	double x;
	double t;
	double told;
	double tout;
	double tend;
	double h;
	double del;
	double absdel;

	// variable sized arrays;
	double* y;
	double* yy;

	// member functions
	char sign(double);

	void test_inputs();
	void setup();
	void first_step();

};

