#include <iostream>
#include "Ode.h"

#define PI 3.14159265358979323846264338328

using namespace std;

template <class T>
int sgn(T a) { return (a > T(0)) - (a < T(0)); }

void f(double x, double* y, double* yp)
{
	double temp = y[0] * y[0] + y[1] * y[1];
	double mu = 1;
	temp *= sqrt(temp);
	yp[0] = y[2];
	yp[1] = y[3];
	yp[2] = -mu / temp * y[0];
	yp[3] = -mu / temp * y[1];
}

int main()
{
	const unsigned int neqn = 4;
	double work[neqn * 21];

	double t = 0;
	double tout = 16 * PI;
	double relerr = 0;
	double abserr = 1e-14;

	double y[neqn] = { 0, 0.4, 2, 0 };



	Ode myDe(f, neqn, y, t, tout, relerr, abserr, work);


	
	myDe.step();
	

	cout << '\n';
	cout << myDe.y[0];
	cout << '\n';
	cout << myDe.y[1];
	cout << '\n';
	cout << myDe.y[2];
	cout << '\n';
	cout << myDe.y[3];
	cout << '\n';
	cout << '\n';
	cout << myDe.abserr;
	cout << '\n';
	cout << myDe.iflag;
	cout << '\n';
	cout << myDe.nostep;

	myDe.tout = 2 * tout;
	myDe.step();

	cout << '\n';
	cout << '\n';
	cout << myDe.y[0];
	cout << '\n';
	cout << myDe.y[1];
	cout << '\n';
	cout << myDe.y[2];
	cout << '\n';
	cout << myDe.y[3];
	cout << '\n';
	cout << '\n';
	cout << myDe.abserr;
	cout << '\n';
	cout << myDe.iflag;
	cout << '\n';
	cout << myDe.nostep;
	cout << '\n';

	return 0;
}
