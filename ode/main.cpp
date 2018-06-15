#include <iostream>
#include <iomanip>
#include <cmath>
#include "Ode.h"
#include "Diffeq.h"

#define PI 3.14159265358979323846264338328

int main()
{
	const int neqn = 42;
	double work[neqn * 22];

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

	Ode myDe(diffeq::seven_body, neqn, y, t, tout, relerr, abserr, work);
	myDe.step();
	
	std::cout << std::fixed;
	std::cout << std::setprecision(5);
	for (int i = 0; i < 7; i++)
	{
		for (int j = 0; j < 6; j++)
		{
			std::cout << std::setw(8) << std::setfill('+') << myDe.y[6 * i + j] << " ";
		}
		std::cout << '\n';
	}
	std::cout << myDe.iflag;
	std::cout << '\n';
	std::cout << myDe.nostep;
	std::cout << '\n';
	std::cin.get();
	return 0;
}
