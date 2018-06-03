#include "Step.h"



Step::Step()
{
	double u = machine();
	twou = 2 * u;
	fouru = 4 * u;
}


Step::~Step()
{
}

void Step::block0()
{
	crash = false;
	test_inputs();
	if (crash) { return; }
	if (start) { initialize(); }

}

char Step::sign(double a)
{
	int b = a > 0 ? 1 : -1;
	return b;
}

double Step::machine()
{
	double halfu = 0.5;
	while (1 + halfu > 1)
	{
		halfu *= 0.5;
	}
	double u = 2 * halfu;
	return u;
}

void Step::test_inputs()
{
	if (abs(h) < fouru * abs(x))
	{
		h = fouru * abs(x) * sign(h);
		crash = true;
		return;
	}
	
	double round = 0.0;
	for (size_t i = 0; i < neqn; i++)
	{
		round += pow(y[i] / wt[i], 2);
	}
	round = twou * sqrt(round);

	if (0.5 * eps < round)
	{
		eps = 2.0 * round * (1 + fouru);
		crash = true;
		return;
	}
}

void Step::initialize()
{
	f(x, y, yp);

}