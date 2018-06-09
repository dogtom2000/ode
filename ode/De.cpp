#include "De.h"
#include <algorithm>

#define sign(a) (a > 0 ? 1 : -1)
#define phi(j, i) phi[j + i * neqn]

using namespace std;

De::De(void(*f)(double, double[], double[]), unsigned char neqn, double* y, double t, double tout, double relerr, double abserr, double* work) : f(f), neqn(neqn), y(y), t(t), tout(tout), relerr(relerr), abserr(abserr)
{
	yp = work + 0 * neqn;
	yy = work + 1 * neqn;
	wt = work + 2 * neqn;
	ypout = work + 3 * neqn;
	yout = work + 4 * neqn;
	phi = work + 5 * neqn;

	iflag = 1;
	machine();

	destep.y = yy;
	destep.y[0] = 5.0;
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
	if ((iflag == 1) || (isnold < 0) || (delsn * del <= 0)) { first_step(); }

	//while (true)
	for (size_t run = 0; run < 3; run++)
	{
		// beyond tout?
		if (abs(x - t) >= absdel) {
			btout();
			return;
		}

		// close to tout and cannot pass?
		if ((isn < 0) || (abs(tout - x) < destep.fouru * abs(x)))
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

void De::machine()
{
	double halfu = 0.5;
	while (1 + halfu > 1)
	{
		halfu *= 0.5;
	}
	 u = 2 * halfu;
	 twou = 2 * u;
	 fouru = 4 * u;
}

void De::test_inputs()
{
	if (neqn > 100) { iflag = 6;  return; }
	if (t == tout) { iflag = 6;  return; }
	if (relerr < 0 || abserr < 0) { iflag = 6;  return; }
	eps = max(relerr, abserr);
	if (eps <= 0) { iflag = 6;  return; }
	if (iflag == 0) { iflag = 6; return; }
	isn = sign(iflag);
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
	tend = 10.0 * del;
	if (isn < 0) { tend = tout; }
	nostep = 0;
	kle4 = 0;
	stiff = false;
	releps = relerr / eps;
	abseps = abserr / eps;
}

void De::first_step()
{
	start = true;
	x = t;
	delsn = sign(del);
	h = max(abs(tout - x), fouru * abs(x)) * sign(tout - x);
	for (size_t l = 0; l < neqn; l++)
	{
		yy[l] = y[l];
	}
}

void De::btout()
{
	interp();
	iflag = 2;
	t = tout;
	isnold = isn;
}

void De::ctout()
{
	h = tout - x;
	f(x, yy, yp);
	for (size_t l = 0; l < neqn; l++)
	{
		y[l] = yy[l] + h * yp[l];
	}
	iflag = 2;
	t = tout;
	told = t;
	isnold = isn;

}

void De::work()
{
	iflag = isn * 4;
	if (stiff) { iflag = isn * 5; }
	for (size_t l = 0; l < neqn; l++)
	{
		y[l] = yy[l];
	}
	t = x;
	told = t;
	isnold = 1;
}

void De::weights()
{
	h = min(abs(h), abs(tend - x)) * sign(h);
	for (size_t l = 0; l < neqn; l++)
	{
		wt[l] = releps * abs(yy[l]) + abseps;
	}
}

void De::tolerances()
{
	iflag = isn * 3;
	relerr = eps * releps;
	abserr = eps * abseps;
	for (size_t l = 0; l < neqn; l++)
	{
		y[l] = yy[l];
	}
	t = x;
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

void De::interp()
{
	double hi = xout - x;
	unsigned char  ki = destep.kold + 1;
	unsigned char  kip1 = ki + 1;

	for (size_t i = 0; i < ki; i++)
	{
		size_t temp1 = i + 1;
		w[i] = 1.0 / temp1;
	}
	double term = 0.0;

	for (size_t j = 1; j < ki; j++)
	{
		int jm1 = j - 1;
		double psijm1 = destep.psi[jm1];
		double gamma = (hi + term) / psijm1;
		double eta = hi / psijm1;
		size_t limit1 = kip1 - (j + 1);
		for (size_t i = 0; i < limit1; i++)
		{
			w[i] = gamma * w[i] - eta * w[i + 1];
		}
		g[j] = w[0];
		rho[j] = gamma * rho[jm1];
		term = psijm1;
	}

	for (size_t l = 0; l < neqn; l++)
	{
		ypout[l] = 0.0;
		yout[l] = 0.0;
	}
	for (size_t j = 0; j < ki; j++)
	{
		int i = kip1 - (j + 1);
		double temp2 = g[i - 1];
		double temp3 = rho[i - 1];
		for (size_t l = 0; l < neqn; l++)
		{
			yout[l] += temp2 * destep.phi(l, i);
			ypout[l] += temp3 * destep.phi(l, i);
		}
	}
	for (size_t l = 0; l < neqn; l++)
	{
		yout[l] = y[l] + h * yout[l];
	}
}