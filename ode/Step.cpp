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
	ifail = 0;
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
	
	double sum = 0.0;
	for (size_t i = 0; i < neqn; i++)
	{
		phi[i * neqn] = yp[i];
		phi[i * neqn + 1] = 0.0;
		sum += pow(yp[i] / wt[i], 2);
	}
	sum = sqrt(sum);
	
	double absh = abs(h);
	if (eps < 16 * sum * h * h)
	{
		absh = 0.25 * sqrt(eps / sum);
	}
	h = fmax(absh, fouru * abs(x)) * sign(h);
	hold = 0.0;
	k = 1;
	start = false;
	phase1 = true;
	nornd = true;

	if (0.5 * eps < 100.0 * round)
	{
		nornd = false;
		for (size_t i = 0; i < neqn; i++)
		{
			phi[i * neqn + 14] = 0.0;
		}
	}
}