#include <iostream>
#include "De.h"
#include "Ode.h"

#define PI 3.14159265358979323846264338328

using namespace std;

template <class T>
int sgn(T a) { return (a > T(0)) - (a < T(0)); }

void f(double x, double* y, double* yp)
{
	double temp = y[0] * y[0] + y[1] * y[1];
	double mu = 398600441800000;
	temp *= sqrt(temp);
	yp[0] = y[2];
	yp[1] = y[3];
	yp[2] = -mu / temp * y[0];
	yp[3] = -mu / temp * y[1];
}

void g(double x, double* y, double* yp)
{
	yp[0] = y[0] * (4 * x * x * x - y[0]) / (x * x * x * x - 1);
}

int main()
{
	const unsigned int neqn = 4;
	// double y[neqn] = { 0, 0.4, 2, 0 };
	double work[neqn * 21];

	double t = 0;
	double tout = 5422.475921;
	tout = 16 * PI;
	double relerr = 0;
	double abserr = 1e-10;

	double y[neqn] = { 0, 6671000, 10850, 0 };
	tout = 1054638.84756323;


	Ode myDe(f, neqn, y, t, tout, relerr, abserr, work);

	//double y2[1] = { 15 };
	//double work2[1 * 21];

	//De myDe(g, 1, y2, 2, 3, 0, 1e-7, work2);


	while (true)
	{
		myDe.step();
		if (myDe.t == tout)
		{
			break;
		}
	}

	
	double aa = { 2 };
	double bb[1] = { 15 };
	double cc[1];

	g(aa, bb, cc);
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
	cout << myDe.integrate.yout[0];
	cout << '\n';
	cout << myDe.integrate.yout[1];
	cout << '\n';
	cout << myDe.integrate.yout[2];
	cout << '\n';
	cout << myDe.integrate.yout[3];
	cout << '\n';
	cout << myDe.integrate.maxk;
	cout << '\n';
	cout << myDe.abserr;
	cout << '\n';
	cout << myDe.relerr;
	return 0;
}
