#include "Step.h"
#include <algorithm>

#define sign(a) (a > 0 ? 1 : -1)
#define aif(a) (a > 0 ? 'p' : a < 0 ? 'n' : '0')
#define p1(a) (a + 1)
#define m1(a) (a - 1)
#define phi(j, i) phi[j + i * neqn]

using namespace std;

Step::Step()
{
}

Step::~Step()
{
}

void Step::take_step()
{
	block0();
}

void Step::block0()
{
	// test if step size and eps are too small
	crash = false;
	test_inputs();
	if (crash) { return; }

	// if this is the first step inialize
	if (start) { initialize(); }

	// number of failed iterations is 0
	ifail = 0;
}

void Step::block1()
{
	// set values offsets
	kp1 = k + 1;
	kp2 = k + 2;
	km1 = k - 1;
	km2 = k - 2;

	// if step size has changed reset step counter
	if (h != hold) { ns = 0; }
	ns = min(ns + 1, kold + 1);
	nsp1 = ns + 1;

	if (k >= ns)
	{
		compute_coefficients();

		if (ns == 1) { initialize_vw(); }
		else { update_vw(); }

		if (k > ns) { compute_g(); }
	}
}

void Step::block2()
{
	if (k > ns) { phi_star(); }

	predict();

	xold = x;
	x += h;
	f(x, p, yp);

	estimate_error();
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
	kold = k;
	hold = h;

	correct();

	f(x, y, yp);

	update_dif();

	double erkp1 = 0.0;
	if ((knew == km1) || (k == 12)) { phase1 = false; }
	if (phase1)
	{
		k = kp1;
		erk = erkp1;
		h *= 2;
	}
	else
	{
		update_h();
	}
}

void Step::test_inputs()
{
	// test if step size is too small, if it is increase it and crash
	if (abs(h) < fouru * abs(x))
	{
		h = fouru * abs(x) * sign(h);
		crash = true;
		return;
	}
	
	// calculate sum to compare to eps in order to control round off error
	double round = 0.0;
	for (size_t l = 0; l < neqn; l++)
	{
		round += pow(y[l] / wt[l], 2);
	}
	round = twou * sqrt(round);

	// test if eps is too small, if it is increase it and crash
	if (0.5 * eps < round)
	{
		eps = 2.0 * round * (1 + fouru);
		crash = true;
		return;
	}
}

void Step::initialize()
{
	// call f to evaluate derivatives
	f(x, y, yp);
	
	// initialize phi
	// calculate sum to compare to h in order to select initial step size
	double sum = 0.0;
	for (size_t l = 0; l < neqn; l++)
	{
		phi(l, 0) = yp[l];
		phi(l, 1) = 0.0;
		sum += pow(yp[l] / wt[l], 2);
	}
	sum = sqrt(sum);
	
	// test if step size is too small, if it is increase it
	absh = abs(h);
	if (eps < 16 * sum * h * h)
	{
		absh = 0.25 * sqrt(eps / sum);
	}
	h = max(absh, fouru * abs(x)) * sign(h);

	// initialize values
	hold = 0.0;
	k = 1;
	kold = 0;
	start = false;
	phase1 = true;
	nornd = true;

	// test if eps is small enough to enable propagated round off control
	if (0.5 * eps < 100.0 * round)
	{
		nornd = false;
		for (size_t l = 0; l < neqn; l++)
		{
			phi(l, 14) = 0.0;
		}
	}
}

void Step::compute_coefficients()
{
	beta[m1(ns)] = 1.0;
	alpha[m1(ns)] = 1.0 / ns;
	double temp1 = h * ns;
	sigma[m1(nsp1)] = 1.0;
	if (k > ns)
	{
		for (size_t i = m1(nsp1); i < k; i++)
		{
			size_t im1 = i - 1;
			double temp2 = psi[im1];
			psi[im1] = temp1;
			beta[i] = beta[im1] * psi[im1] / temp2;
			temp1 = temp2 + h;
			alpha[i] = h / temp1;
			sigma[i + 1] = p1(i) * alpha[i] * sigma[i];
		}
	}
	psi[m1(k)] = temp1;
}

void Step::initialize_vw()
{
	for (size_t iq = 0; iq < k; iq++)
	{
		double temp3 = p1(iq) * (p1(iq) + 1);
		v[iq] = 1.0 / temp3;
		w[iq] = v[iq];
	}
}

void Step::update_vw()
{
	if (k > kold)
	{
		double temp4 = k * kp1;
		v[m1(k)] = 1.0 / temp4;
		nsm2 = ns - 2;
		if (ns >= 3)
		{
			for (size_t j = 0; j < nsm2; j++)
			{
				size_t i = k - p1(j);
				v[m1(i)] -= alpha[j + 1] * v[m1(i) + 1];
			}
		}
	}

	size_t limit1 = kp1 - ns;
	double temp5 = alpha[m1(ns)];
	for (size_t iq = 0; iq < limit1; iq++)
	{
		v[iq] -= temp5 * v[iq + 1];
		w[iq] = v[iq];
	}
	g[m1(nsp1)] = w[0];
}

void Step::compute_g()
{
	nsp2 = ns + 2;
	for (size_t i = m1(nsp2); i < kp1; i++)
	{
		size_t limit2 = kp2 - p1(i);
		double temp6 = alpha[i - 1];
		for (size_t iq = 0; iq < limit2; iq++)
		{
			w[iq] -= temp6 * w[iq + 1];
		}
		g[i] = w[0];
	}
}

void Step::phi_star() 
{
	for (size_t i = m1(nsp1); i < k; i++)
	{
		double temp1 = beta[i];
		for (size_t l = 0; l < neqn; l++)
		{
			phi(l, i) *= temp1;
		}
	}
}

void Step::predict()
{
	for (size_t l = 0; l < neqn; l++)
	{
		phi(l, m1(kp2)) = phi(l, m1(kp1));
		phi(l, m1(kp1)) = 0.0;
		p[l] = 0.0;
	}

	for (size_t j = 0; j < k; j++)
	{
		size_t i = m1(kp1 - p1(j));
		size_t ip1 = i + 1;
		double temp2 = g[i];
		for (size_t l = 0; l < neqn; l++)
		{
			p[l] += temp2 * phi(l, i);
			phi(l, i) += phi(l, ip1);
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
			double tau = h * p[l] - phi(l, 14);
			p[l] = y[l] + tau;
			phi(l, 15) = (p[l] - y[l]) - tau;
		}
	}
}

void Step::estimate_error()
{
	absh = abs(h);
	erkm2 = 0.0;
	erkm1 = 0.0;
	erk = 0.0;

	for (size_t l = 0; l < neqn; l++)
	{
		double temp3 = 1.0 / wt[l];
		double temp4 = yp[l] - phi(l, 0);
		switch (aif(km2))
		{
		case 'p':
			erkm2 += pow((phi(l, m1(km1)) + temp4) * temp3, 2);
		case '0':
			erkm1 += pow((phi(l, m1(k)) + temp4) * temp3, 2);
		case 'n':
			erk += pow(temp4 * temp3, 2);
		}
	}

	switch (aif(km2))
	{
	case 'p':
		erkm2 = absh * sigma[m1(km1)] * gstar[m1(km2)] * sqrt(erkm2);
	case '0':
		erkm1 = absh * sigma[m1(k)] * gstar[m1(km1)] * sqrt(erkm1);
	case 'n':
		break;
	}

	double temp5 = absh * sqrt(erk);
	double err = temp5 * (g[m1(k)] - g[m1(kp1)]);
	erk = temp5 * sigma[m1(kp1)] * gstar[m1(k)];
	knew = k;

	if (km2 > 0)
	{
		if (max(erkm1, erkm2) < erk || (erkm1 <= 0.5 * erk))
		{
			knew = km1;
		}
	}
}

void Step::restore()
{
	for (size_t i = 0; i < k; i++)
	{
		double temp1 = 1.0 / beta[i];
		size_t ip1 = i + 1;
		for (size_t l = 0; l < neqn; l++)
		{
			phi(l, i) = temp1 * (phi(l, i) - phi(l, ip1));
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

void Step::correct()
{
	double temp1 = h * g[m1(k)];
	if (nornd)
	{
		for (size_t l = 0; l < neqn; l++)
		{
			y[l] = p[l] + temp1 * (yp[l] - phi(l, 0));
		}
	}
	else
	{
		for (size_t l = 0; l < neqn; l++)
		{
			double rho = temp1 * (yp[l] - phi(l, 0)) - phi(l, 15);
			y[l] = p[l] + rho;
			phi(l, 14) = (y[l] - p[l]) - rho;
		}
	}
}

void Step::update_dif()
{
	for (size_t l = 0; l < neqn; l++)
	{
		phi(l, m1(kp1)) = yp[l] - phi(l, 0);
		phi(l, m1(kp2)) = phi(l, m1(kp1)) - phi(l, m1(kp2));
	}
	for (size_t i = 0; i < k; i++)
	{
		for (size_t l = 0; l < neqn; l++)
		{
			phi(l, i) += phi(l, m1(kp1));
		}
	}
}

void Step::update_h()
{
	if (knew == km1)
	{
		k = km1;
		erk = erkm1;
	}
	else if (k < ns)
	{
		for (size_t l = 0; l < neqn; l++)
		{
			erkp1 += pow(phi(l, m1(kp2)) / wt[l], 2);
		}
		erkp1 = absh * gstar[m1(kp1)] * sqrt(erkp1);
		if (k == 1)
		{
			if (erkp1 < 0.5 * erk)
			{
				k = kp1;
				erk = erkp1;
			}
		}
		else
		{
			if (erkm1 <= min(erk, erkp1))
			{
				k = km1;
				erk = erkm1;
			}
			else if ((erkp1 >= erk) || (k == 12))
			{
			}
			else
			{
				k = kp1;
				erk = erkp1;
			}
		}
	}



	// if knew == km1 goto 455

	// if kp1 > ns goto 460
	
	// else do loop

	// if k > 1 goto 445

	// if erkp1 >= 0.5 * erk goto 460

	// else goto 450

	// 445 goto 450, 455, 460

	// 450 goto 460

	// 455 goto 460




	// 460
	double hnew = h + h;
	if (0.5 * eps >= erk * two[k])
	{
		h = hnew;
		return;
	}
	hnew = h;
	if (0.5 * eps >= erk)
	{
		h = hnew;
		return;
	}
	size_t temp2 = k + 1;
	double r = pow((0.5 * eps / erk), 1.0 / temp2);
	hnew = absh * max(0.5, min(0.9, r));
	hnew = max(hnew, fouru * abs(x)) * sign(h);
	h = hnew;
	return;
}