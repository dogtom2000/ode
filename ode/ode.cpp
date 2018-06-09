#include <iostream>
#include "De.h"
#define sign(a) (a > 0 ? 1 : -1)

using namespace std;


void f(double x, double* y, double* yp)
{
	yp[0] = y[2];
	yp[1] = y[3];
	yp[2] = -3.986e14 / abs(y[0] * y[0]) * sign(y[0]);
	yp[3] = -3.986e14 / abs(y[1] * y[1]) * sign(y[0]);
}

int main()
{
	const unsigned char neqn = 4;
	double y[neqn] = { 0, 6671000, 7800, 0 };
	double work[neqn * 21];

	double t = 0;
	double tout = 60 * 90;
	double relerr = 1e-5;
	double abserr = 1e-8;

	De myDe(f, neqn, y, t, tout, relerr, abserr, 1, work);
	myDe.step();

	cout << myDe.wt[0];
	cout << '\n';
	cout << myDe.yy[1];
	cout << '\n';
	cout << myDe.yy[2];
	cout << '\n';
	cout << myDe.yy[3];
	cout << '\n';

	return 0;
}