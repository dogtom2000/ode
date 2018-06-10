#include "De.h"

#define phi(j, i) phi[j + i * neqn]

using namespace std;

De::De(void(*f)(double, double*, double*), unsigned char neqn, double* y, double t, double tout, double relerr, double abserr, double* work)
	: neqn(neqn), y(y), t(t), tout(tout), relerr(relerr), abserr(abserr)
{
	iflag = 1;
	
	destep.y		= work + 0 * neqn;
	destep.yp		= work + 1 * neqn;
	destep.yout		= work + 2 * neqn;
	destep.ypout	= work + 3 * neqn;
	destep.p		= work + 4 * neqn;
	destep.wt		= work + 5 * neqn;
	destep.phi		= work + 6 * neqn;
	destep.phi1 = work + 6 * neqn + 1;
	destep.phi2 = work + 6 * neqn + 2;
	destep.phi3 = work + 6 * neqn + 3;
	destep.phi4 = work + 6 * neqn + 4;
	destep.phi5 = work + 6 * neqn + 5;

	destep.f = f;

	double u = machine();
	destep.twou = 2 * u;
	destep.fouru = 4 * u;
	destep.neqn = neqn;
}


De::~De()
{
}

void De::step()
{
	// test for valid inputs
	test_inputs();
	if (iflag == 6) { return; }

	// set up integration parameters
	setup();

	// determine if work quantities need set
	if ((iflag == 1) || (isnold < 0) || (delsgn * del <= 0)) { first_step(); }

	while (true)
	{
		// beyond tout?
		if (abs(destep.x - t) >= absdel) {
			btout();
			return;
		}

		// close to tout and cannot pass?
		if ((isgn < 0) || (abs(tout - destep.x) < destep.fouru * abs(destep.x)))
		{
			ctout();
			return;
		}

		// test for too much work
		if (nostep > maxnum)
		{
			work();
			return;
		}

		weights();
		destep.take_step();
		if (destep.crash)
		{
			tolerances();
			return;
		}
		else
		{
			increment();
			
		}
	}
}

double De::machine()
{
	double halfu = 0.5;
	while (1 + halfu > 1)
	{
		halfu *= 0.5;
	}
	double u = 2 * halfu;
	return u;
}

void De::test_inputs()
{
	if (neqn > 100) { iflag = 6;  return; }
	if (t == tout) { iflag = 6;  return; }
	if (relerr < 0 || abserr < 0) { iflag = 6;  return; }
	destep.eps = max(relerr, abserr);
	if (destep.eps <= 0) { iflag = 6;  return; }
	if (iflag == 0) { iflag = 6; return; }
	isgn = sgn(iflag);
	iflag = abs(iflag);
	if (iflag == 1)
	{
		return;
	}
	else
	{
		if (t != told) { iflag = 6; return; }
	}
	if (iflag < 2 || iflag > 5) { iflag = 6; return; }
}

void De::setup()
{
	del = tout - t;
	absdel = abs(del);
	tend = t + 10.0 * del;
	if (isgn < 0) { tend = tout; }
	nostep = 0;
	kle4 = 0;
	stiff = false;
	releps = relerr / destep.eps;
	abseps = abserr / destep.eps;
}

void De::first_step()
{
	destep.start = true;
	destep.x = t;
	delsgn = sgn(del);
	destep.h = max(abs(tout - destep.x), destep.fouru * abs(destep.x)) * sgn(tout - destep.x);
	for (size_t l = 0; l < neqn; l++)
	{
		destep.y[l] = y[l];
	}
}

void De::btout()
{
	destep.xout = tout;
	destep.interp();
	iflag = 2;
	t = tout;
	isnold = isgn;
}

void De::ctout()
{
	destep.h = tout - destep.x;
	destep.extrap();
	iflag = 2;
	t = tout;
	told = t;
	isnold = isgn;
}

void De::work()
{
	iflag = isgn * 4;
	if (stiff) { iflag = isgn * 5; }
	for (size_t l = 0; l < neqn; l++)
	{
		y[l] = destep.y[l];
	}
	t = destep.x;
	told = t;
	isnold = 1;
}

void De::weights()
{
	destep.h = min(abs(destep.h), abs(tend - destep.x)) * sgn(destep.h);
	for (size_t l = 0; l < neqn; l++)
	{
		destep.wt[l] = releps * abs(destep.y[l]) + abseps;
	}
}

void De::tolerances()
{
	iflag = isgn * 3;
	relerr = destep.eps * releps;
	abserr = destep.eps * abseps;
	for (size_t l = 0; l < neqn; l++)
	{
		y[l] = destep.y[l];
	}
	t = destep.x;
	told = t;
	isnold = 1;
}

void De::increment()
{
	nostep++;
	kle4++;
	if (destep.kold > 4) { kle4 = 0; }
	if (kle4 > 50) { stiff = true; }

}