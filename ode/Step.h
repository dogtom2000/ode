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
	bool phase1;
	bool nornd;


	// variables
	unsigned char neqn;
	unsigned char k;
	unsigned char ifail;

	double twou;
	double fouru;
	double eps;
	double round;

	double h;
	double hold;
	double x;

	

	// variable sized arrays;
	double* y;
	double* yp;
	double* wt;
	double* phi;

	// fixed size arrays

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

