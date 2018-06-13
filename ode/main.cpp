#include <iostream>
#include <cmath>
#include "Ode.h"

#define PI 3.14159265358979323846264338328

template <class T>
int sgn(T a) { return (a > T(0)) - (a < T(0)); }

void f(double x, double* y, double* yp)
{
	double temp = y[0] * y[0] + y[1] * y[1];
	double mu = 1;
	temp *= std::sqrt(temp);
	yp[0] = y[2];
	yp[1] = y[3];
	yp[2] = -mu / temp * y[0];
	yp[3] = -mu / temp * y[1];
}

void g(double x, double* y, double* yp)
{
	const int nbody = 7;
	const double mu[nbody] = { 1.0, 0.00000244784, 0.000003003489663, 0.000000036958626, 0.0000003227155, 0.000954594263, 0.00028581485 };

	for (size_t i = 0; i < nbody; i++)
	{
		int ii = 6 * i;
		for (size_t j = 0; j < 3; j++)
		{
			yp[ii + j] = y[ii + j + 3];
		}

		yp[ii + 3] = 0.0;
		yp[ii + 4] = 0.0;
		yp[ii + 5] = 0.0;
		for (size_t j = 0; j < nbody; j++)
		{
			if (i != j)
			{
			int jj = 6 * j;
			double rx = y[ii + 0] - y[jj + 0];
			double ry = y[ii + 1] - y[jj + 1];
			double rz = y[ii + 2] - y[jj + 2];
			double r = std::sqrt(rx * rx + ry * ry + rz * rz);

			yp[ii + 3] -= mu[j] * rx / std::pow(r, 3);
			yp[ii + 4] -= mu[j] * ry / std::pow(r, 3);
			yp[ii + 5] -= mu[j] * rz / std::pow(r, 3);
			}
		}
	}
}

int main()
{
	//const unsigned int neqn = 4;
	//double work[neqn * 21];

	//double t = 0;
	//double tout = 16 * PI;
	//double relerr = 0;
	//double abserr = 1e-14;

	//double y[neqn] = { 0, 0.4, 2, 0 };

	const unsigned int neqn = 42;
	double work[neqn * 21];

	double t = 0;
	double tout = 2.0 * PI;
	double relerr = 0;
	double abserr = 2e-8;

	double y[neqn]{	0,			0, 0, 0, 0,			0,
					0.723,		0, 0, 0, 1.1760637,	0,
					1,			0, 0, 0, 1,			0,
					1.00257,	0, 0, 0, 1,			0,
					1.524,		0, 0, 0, 0.810042,	0,
					5.203,		0, 0, 0, 0.4384026,	0,
					9.537,		0, 0, 0, 0.3238129,	0 };








	Ode myDe(g, neqn, y, t, tout, relerr, abserr, work);


	
	myDe.step();
	
	
	//std::cout << '\n';
	//std::cout << myDe.y[0];
	//std::cout << '\n';
	//std::cout << myDe.y[1];
	//std::cout << '\n';
	//std::cout << myDe.y[2];
	//std::cout << '\n';
	//std::cout << myDe.y[3];
	//std::cout << '\n';
	//std::cout << '\n';
	std::cout << myDe.abserr;
	std::cout << '\n';
	std::cout << myDe.iflag;
	std::cout << '\n';
	std::cout << myDe.nostep;


	return 0;
}
