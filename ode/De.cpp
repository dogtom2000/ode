#include "De.h"



De::De()
{
}


De::~De()
{
}

char De::sign(double a)
{
	int b = a > 0 ? 1 : -1;
	return b;
}

void De::test_inputs()
{
	if (neqn < 0 || neqn > 20) { iflag = 6;  return; }
	if (t == tout) { iflag = 6;  return; }
	if (relerr < 0 || abserr < 0) { iflag = 6;  return; }
	if (eps <= 0) { iflag = 6;  return; }
	if (iflag == 0) { iflag = 6; return; }
	isn = sign(iflag);
	iflag = abs(iflag);
	if (iflag != 1)
	{
		if (t != told) { iflag = 6; return; }
	}
	if (!(iflag >= 2 && iflag <=5)) { iflag = 6; return; }
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
	h = fmax(abs(t - x), destep.fouru * abs(x)) * sign(tout - x);
	for (size_t l = 0; l < neqn; l++)
	{
		yy[l] = y[l];
	}
}