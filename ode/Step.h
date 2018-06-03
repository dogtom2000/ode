#pragma once
#include <cmath>
class Step
{
public:
	// constants

	// constructor and destructor
	Step();
	~Step();

	// logical
	bool crash;
	bool start;


	// variables
	unsigned char neqn = 1;

	double twou;
	double fouru;
	double eps;

	double h;
	double x;

	

	// arrays
	double y[1];
	double yp[1];
	double wt[1];

	// member functions
	void(*f)(double, double[], double[]);

	char sign(double);
	double machine();
	
	void block0();
	//void block1();
	//void block2();
	//void block3();
	//void block4();

	void test_inputs();
	void initialize();
	//void compute_coefficients();
	//void initialize_vw();
	//void update_vw();
	//void compute_g();

};

