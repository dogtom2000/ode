#include "Ode.h"
#include <iostream>

Ode::Ode(void(*f)(double, double*, double*), unsigned int neqn, double* y, double t, double tout, double relerr, double abserr, double* work)
	: neqn(neqn), y(y), t(t), tout(tout), relerr(relerr), abserr(abserr)
{
	iflag = 1;

	integrate.y = work + 0 * neqn;
	integrate.yp = work + 1 * neqn;
	integrate.yout = work + 2 * neqn;
	integrate.ypout = work + 3 * neqn;
	integrate.p = work + 4 * neqn;
	integrate.wt = work + 5 * neqn;
	integrate.phi = work + 6 * neqn;
	integrate.phi1 = work + 6 * neqn + 1;
	integrate.phi2 = work + 6 * neqn + 2;
	integrate.phi3 = work + 6 * neqn + 3;
	integrate.phi4 = work + 6 * neqn + 4;
	integrate.phi5 = work + 6 * neqn + 5;

	integrate.neqn = neqn;
	integrate.f = f;

	integrate.maxk = 1;

	double halfu = 0.5;
	while (1 + halfu > 1)
	{
		halfu *= 0.5;
	}
	integrate.twou = 4 * halfu;
	integrate.fouru = 8 * halfu;
}

Ode::~Ode()
{
}

void Ode::step()
{
	if (!test_inputs())
	{
		return;
	}

	setup();

	if ((iflag == 1) || (isgnold < 0) || (delsgn * del <= 0))
	{ 
		first_step();
	}

	while (nostep < maxnum)
	{
		if (abs(integrate.x - t) >= absdel)
		{
			end_interp();
			return;
		}

		if ((isgn < 0) && (abs(tout - integrate.x) < integrate.fouru * abs(integrate.x)))
		{
			end_extrap();
			return;
		}

		set_weights();

		integrate.take_step();

		std::cout << integrate.k;
		std::cout << ' ';

		if (integrate.crash)
		{			
			end_tol();
			return;
		}
		else
		{
			increment();
		}
	}

	end_work();
	return;
}

// incomplete
bool Ode::test_inputs()
{
	integrate.eps = std::max(relerr, abserr);
	isgn = sgn(iflag);
	iflag = abs(iflag);

	return true;
}

void Ode::setup()
{
	del = tout - t;
	absdel = abs(del);
	tend = t + 10.0 * del;

	if (isgn < 0) { tend = tout; }
	nostep = 0;
	kle4 = 0;
	stiff = false;
	releps = relerr / integrate.eps;
	abseps = abserr / integrate.eps;
}

void Ode::first_step()
{
	integrate.start = true;
	integrate.x = t;
	delsgn = sgn(del);
	integrate.h = sgn(tout - integrate.x) * std::max(abs(tout - integrate.x), integrate.fouru * abs(integrate.x));
	for (size_t l = 0; l < neqn; l++)
	{
		integrate.y[l] = y[l];
	}
}

void Ode::set_weights()
{
	integrate.h = sgn(integrate.h) * std::min(abs(integrate.h), abs(tend - integrate.x));
	for (size_t l = 0; l < neqn; l++)
	{
		integrate.wt[l] = releps * abs(integrate.y[l]) + abseps;
	}
}

void Ode::increment()
{
	nostep++;
	kle4++;
	if (integrate.kold > 4) { kle4 = 0; }
	if (kle4 > 50) { stiff = true; }
}

void Ode::end_interp()
{
	integrate.xout = tout;
	integrate.interp();
	iflag = 2;
	t = tout;
	told = t;
	isgnold = isgn;
}

void Ode::end_extrap()
{
	integrate.h = tout - integrate.x;
	integrate.extrap();
	iflag = 2;
	t = tout;
	told = t;
	isgnold = isgn;
}

void Ode::end_work()
{
	iflag = isgn * 4;
	if (stiff) { iflag = isgn * 5; }
	for (size_t l = 0; l < neqn; l++)
	{
		y[l] = integrate.y[l];
	}
	t = integrate.x;
	told = t;
	isgnold = 1;
}

void Ode::end_tol()
{
	iflag = isgn * 3;
	relerr = integrate.eps * releps;
	abserr = integrate.eps * abseps;
	for (size_t l = 0; l < neqn; l++)
	{
		y[l] = integrate.y[l];
	}
	t = integrate.x;
	told = t;
	isgnold = 1;
}

