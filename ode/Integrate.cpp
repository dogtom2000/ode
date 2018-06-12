#include "Integrate.h"
#include <algorithm>

#define phi(j, i) phi[j + i * neqn]


Integrate::Integrate()
{
	gi[0] = 1;
	rho[0] = 1;
}


Integrate::~Integrate()
{
}


void Integrate::take_step()
{
	block0();
	if (crash) { return; }

	while (true)
	{
		block1();
		block2();
		if (step_fail)
		{
			block3();
			if (crash) { return; }
		}
		else
		{
			block4();
			return;
		}
	}
}

void Integrate::interp()
{
	double hi = xout - x;
	unsigned int ki = kold + 1;
	unsigned int kip1 = ki + 1;

	for (size_t i = 0; i < ki; i++)
	{
		unsigned int temp1 = i + 1;
		wi[i] = 1.0 / temp1;
	}

	double term = 0.0;
	
	for (size_t j = 1; j < ki; j++)
	{
		unsigned int jm1 = j - 1;
		double psijm1 = psi[jm1];
		double gamma = (hi + term) / psijm1;
		double eta = hi / psijm1;
		unsigned int limit1 = kip1 - j - 1;
		for (size_t i = 0; i < limit1; i++)
		{
			wi[i] = gamma * wi[i] - eta * wi[i + 1];
		}
		gi[j] = wi[0];
		rho[j] = gamma * rho[jm1];
		term = psijm1;
	}

	for (size_t l = 0; l < neqn; l++)
	{
		yout[l] = 0.0;
		ypout[l] = 0.0;
	}

	for (size_t j = 0; j < ki; j++)
	{
		unsigned int i = kip1 - j - 2;
		double temp2 = gi[i];
		double temp3 = rho[i];
		for (size_t l = 0; l < neqn; l++)
		{
			yout[l] += temp2 * phi(l, i);
			ypout[l] += temp3 * phi(l, i);
		}
	}

	for (size_t l = 0; l < neqn; l++)
	{
		yout[l] = y[l] + hi * yout[l];
	}
}

void Integrate::extrap()
{

}

void Integrate::block0()
{

}

void Integrate::block1()
{

}

void Integrate::block2()
{

}

void Integrate::block3()
{

}

void Integrate::block4()
{

}

// block 0 functions
void Integrate::test_inputs()
{

}

void Integrate::setup()
{

}

// block 1 functions
void Integrate::compute_coefficients()
{

}

void Integrate::initialize_vw()
{

}

void Integrate::update_vw()
{

}

void Integrate::compute_g()
{

}

// block 2 functions
void Integrate::phi_star()
{

}

void Integrate::predict()
{

}

void Integrate::estimate_error()
{

}

// block 3 functions
void Integrate::restore()
{

}

void Integrate::order_one()
{

}

// block 4 functions
void Integrate::correct()
{

}

void Integrate::update_dif()
{

}

void Integrate::update_h()
{

}