#pragma once
#include <algorithm>
#include <cmath>
#include "Step.h"

class De
{
public:
	// constants
	const short maxnum{ 500 };

	// constructor and destructor
	De(void(*)(double, double*, double*), unsigned char, double*, double, double, double, double, double*);
	~De();

	// logical
	bool stiff;

	// variables
	char iflag;
	char isgn;
	char isnold;
	char delsgn;
	
	unsigned char neqn;
	unsigned char nostep;
	unsigned char kle4;

	double relerr;
	double abserr;
	double releps;
	double abseps;
	double t;
	double told;
	double tout;
	double tend;
	double del;
	double absdel;

	// variable sized arrays;
	double* y;
	double* yout;

	// step object
	Step destep;

	// member functions
	double machine();

	void step();

	void test_inputs();
	void setup();
	void first_step();
	void btout();
	void ctout();
	void work();
	void weights();
	void tolerances();
	void increment();

	// inline functions
	template <class T>
	int sgn(T a) { return (a > T(0)) - (a < T(0)); }
};

