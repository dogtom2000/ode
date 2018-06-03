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

	predict1();

	xold = x;
	x += h;
	f(x, p, yp);

	double err = estimate_error();
	step_success = err <= eps ? true : false;
}

void Step::block3()
{
	phase1 = false;
	x = xold;

	restore();

	ifail++;

	order_one();
}

void Step::block4()
{

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
	for (size_t l = 0; l < neqn; l++)
	{
		round += pow(y[l] / wt[l], 2);
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
	for (size_t l = 0; l < neqn; l++)
	{
		phi[0 * neqn + l] = yp[l];
		phi[1 * neqn + l] = 0.0;
		sum += pow(yp[l] / wt[l], 2);
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
		for (size_t l = 0; l < neqn; l++)
		{
			phi[14 * neqn + l] = 0.0;
		}
	}
}

void Step::compute_coefficients()
{
	beta[ns - 1] = 1.0;
	alpha[ns - 1] = 1.0 / ns;
	double temp1 = h * ns;
	sigma[ns] = 1.0;
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
	for (size_t iq = 0; iq < k; iq++)
	{
		double temp3 = (iq + 1) * (iq + 2);
		v[iq] = 1.0 / temp3;
		w[iq] = v[iq];
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
	for (size_t iq = 0; iq < limit1; iq++)
	{
		v[iq] -= temp5 * v[iq + 1];
		w[iq] = v[iq];
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
			for (size_t iq = 0; iq < limit2; iq++)
			{
				w[iq] -= temp6 * w[iq + 1];
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
		for (size_t l = 0; l < neqn; l++)
		{
			phi[i * neqn + l] *= temp1;
		}
	}
}

void Step::predict1()
{
	for (size_t l = 0; l < neqn; l++)
	{
		phi[kp1 * neqn + l] = phi[kp1 * neqn + l];
		phi[k * neqn + l] = 0.0;
		p[l] = 0.0;
	}
	for (size_t j = 0; j < k; j++)
	{
		size_t i = k - j - 1;
		double temp2 = g[i];
		for (size_t l = 0; l < neqn; l++)
		{
			p[l] += temp2 * phi[i * neqn + l];
			phi[i * neqn + l] += phi[(i + 1) * neqn + l];
		}
	}
	if (nornd)
	{
		for (size_t l = 0; l < neqn; l++)
		{
			p[l] = y[l] + h * p[l];
		}
	}
	else
	{
		for (size_t l = 0; l < neqn; l++)
		{
			double tau = h * p[l] - phi[14 * neqn + l];
			p[l] = y[l] + tau;
			phi[15 * neqn + l] = (p[l] - y[l]) - tau;
		}
	}
}

double Step::estimate_error()
{
	double absh = abs(h);
	double erkm2 = 0.0;
	double erkm1 = 0.0;
	erk = 0.0;
	for (size_t l = 0; l < neqn; l++)
	{
		double temp3 = 1.0 / wt[l];
		double temp4 = yp[l] - phi[0 * neqn + l];
		if (km2 > 0)
		{
			erkm2 += pow(((phi[km2 * neqn + l] + temp4) * temp3), 2);
		} 
		if (km2 >= 0)
		{
			erkm1 += pow(((phi[km1 * neqn + l] + temp4) * temp3), 2);
		}
		erk += pow(temp3 * temp4, 2);
	}
	if (km2 > 0)
	{
		erkm2 = absh * sigma[km2] * gstar[k - 3] * sqrt(erkm2);
	}
	if (km2 >= 0)
	{
		erkm1 = absh * sigma[km1] * gstar[km2] * sqrt(erkm1);
	}
	double temp5 = absh * sqrt(erk);
	double err = temp5 * (g[km1] - g[k]);
	erk = temp5 * sigma[k] * gstar[km1];
	knew = k;
	if (km2 > 0)
	{
		if (fmax(erkm1, erkm2) < erk)
		{
			knew = km1;
		}
	}
	else if (km2 >= 0)
	{
		if (erkm1 <= 0.5 * erk)
		{
			knew = km1;
		}
	}

	return err;
}

void Step::restore()
{
	for (size_t i = 0; i < k; i++)
	{
		double temp1 = 1.0 / beta[i];
		size_t ip1 = i + 1;
		for (size_t l = 0; l < neqn; l++)
		{
			phi[i * neqn + l] = temp1 * (phi[i * neqn + l] - phi[ip1 * neqn + l]);
		}
		if (k > 1)
		{
			for (size_t i = 1; i < k; i++)
			{
				psi[i - 1] = psi[i] - h;
			}
		}
	}
}

void Step::order_one()
{
	double temp2 = 0.5;
	if (ifail > 3)
	{
		if (0.5 * eps < 0.25 * erk)
		{
			temp2 = sqrt(0.5 * eps / erk);
		}
	}
	if (ifail >= 3)
	{
		knew = 1;
	}
	h *= temp2;
	k = knew;
	if (abs(h) < fouru * abs(x))
	{
		crash = true;
		h = fouru * abs(x) * sign(h);
		eps *= 2;
	}
}