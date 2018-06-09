#include <iostream>
#include "De.h"

using namespace std;

template <class T>
int sgn(T a) { return (a > T(0)) - (a < T(0)); }

void f(double x, double* y, double* yp)
{
	double temp = y[0] * y[0] + y[1] * y[1];
	temp *= sqrt(temp);
	yp[0] = y[2];
	yp[1] = y[3];
	yp[2] = -3.986e14 / temp * y[0];
	yp[3] = -3.986e14 / temp * y[1];
}

int main()
{
	const unsigned char neqn = 4;
	double y[neqn] = { 0, 6671000, 7729.887564, 0 };
	double work[neqn * 21];

	double t = 0;
	double tout = 5422.475921;
	double relerr = 1e-14;
	double abserr = 1e-20;

	De myDe(f, neqn, y, t, tout, relerr, abserr, work);
	myDe.step();

	cout << myDe.destep.yout[0];
	cout << '\n';
	cout << myDe.destep.yout[1];
	cout << '\n';
	cout << myDe.destep.yout[2];
	cout << '\n';
	cout << myDe.destep.yout[3];
	cout << '\n';

	return 0;
}