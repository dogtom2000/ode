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
	bool step_success;


	// variables
	char neqn;
	char ns;
	char k;
	char knew;
	char kold;
	char kp1;
	char kp2;
	char km1;
	char km2;
	char ifail;

	double twou;
	double fouru;
	double eps;
	double round;
	double erk;
	double erkm1;
	double erkm2;

	double h;
	double hold;
	double absh;
	double x;
	double xold;

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
	void calculate_two();
	
	void block0();
	void block1();
	void block2();
	void block3();
	void block4();

	// block 0 functions
	void test_inputs();
	void initialize();

	// block 1 functions
	void compute_coefficients();
	void initialize_vw();
	void update_vw();
	void compute_g();

	// block 2 functions
	void phi_star();
	void predict();
	double estimate_error();

	// block 3 functions
	void restore();
	void order_one();

	// block 4 functions
	void correct();
	void update_dif();
	void update_h();
};