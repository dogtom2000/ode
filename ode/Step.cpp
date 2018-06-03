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

void Step::block1()
{
	kp1 = k + 1;
	kp2 = k + 2;
	km1 = k - 1;
	km2 = k - 2;

	if (h != hold) { ns = 0; }
	ns = fmin(ns + 1, k + 1);
	if (k >= ns)
	{
		compute_coefficients();
		if (ns == 1)
		{
			initialize_vw();
		}
		else
		{
			update_vw();
		}
		if (k > ns)
		{
			compute_g();
		}
	}
}

void Step::block2()
{
	if (k > ns)
	{
		phi_star();
	}
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
		phi[i * neqn + 0] = yp[i];
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

void Step::compute_coefficients()
{
	beta[ns - 1] = 1.0;
	alpha[ns - 1] = 1.0 / ns;
	double temp1 = h * ns;
	sigma[ns] = 1;
	if (k >= ns)
	{
		for (size_t i = ns; i < k; i++)
		{
			size_t im1 = i - 1;
			double temp2 = psi[im1];
			psi[im1] = temp1;
			beta[i] = beta[im1] * psi[im1] / temp2;
			temp1 = temp2 + h;
			alpha[i] = h / temp1;
			sigma[i + 1] = (i + 1) * alpha[i] * sigma[i];

		}
	}
	psi[k - 1] = temp1;
}

void Step::initialize_vw()
{
	for (size_t i = 0; i < k; i++)
	{
		double temp3 = (i + 1) * (i + 2);
		v[i] = 1.0 / temp3;
		w[i] = v[i];
	}
}

void Step::update_vw()
{
	if (k > kold)
	{
		double temp4 = k * kp1;
		v[k - 1] = 1.0 / temp4;
		if (ns >= 3)
		{
			for (size_t j = 0; j < ns - 2; j++)
			{
				size_t i = k - j - 2;
				v[i] -= alpha[j + 1] * v[i + 1];
			}
		}
	}

	size_t limit1 = kp1 - ns;
	double temp5 = alpha[ns - 1];
	for (size_t i = 0; i < limit1; i++)
	{
		v[i] -= temp5 * v[i + 1];
		w[i] = v[i];
	}
	g[ns] = w[0];
}

void Step::compute_g()
{
	if (k > ns)
	{
		for (size_t i = ns + 1; i < kp1; i++)
		{
			size_t limit2 = kp1 - i;
			double temp6 = alpha[i - 1];
			for (size_t j = 0; j < limit2; j++)
			{
				w[j] -= temp6 * w[j + 1];
			}
			g[i] = w[0];
		}
	}
}

void Step::phi_star() 
{
	for (size_t i = ns; i < k; i++)
	{
		double temp1 = beta[i];
		for (size_t j = 0; j < neqn; j++)
		{
			phi[i * neqn + j] *= temp1;
		}
	}
}

void Step::predict1()
{
	for (size_t i = 0; i < neqn; i++)
	{
		phi[i * neqn + kp1] = phi[i * neqn + k];
		phi[i * neqn + k] = 0.0;
		p[i] = 0.0;
	}

}

void Step::estimate_error()
{

}

