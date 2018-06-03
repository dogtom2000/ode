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
	size_t neqn;
	size_t ns;
	size_t k;
	size_t kold;
	size_t kp1;
	size_t kp2;
	size_t km1;
	size_t km2;
	size_t ifail;

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
	double* p;
	double* phi;

	// fixed size arrays
	double alpha[12];
	double beta[12];
	double psi[12];
	double sigma[13];

	double v[12];
	double w[12];
	double g[13];
	double gstar[13];
	double two[13];

	// member functions
	void(*f)(double, double[], double[]);

	char sign(double);
	double machine();
	
	void block0();
	void block1();
	void block2();
	//void block3();
	//void block4();

	void test_inputs();
	void initialize();
	void compute_coefficients();
	void initialize_vw();
	void update_vw();
	void compute_g();

	void phi_star();
	void predict1();
	void estimate_error();

};

