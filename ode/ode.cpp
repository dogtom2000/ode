#include <iostream>
#include "Step.h"


double test(double a[][3]) {

	double b = a[1][2];
	return b;
}



int main()
{
	double arr1[2][3] =
	{
		{0, 1, 2},
		{4, 5, 6 }
	};

	double out = test(arr1);

	double arr2[5] = { 0, 1, 2, 3, 4 };

	Step myStep;
	myStep.y = arr2;

	out = myStep.y[3];

	int ns = 4;
	int nsp1 = 5;
	int nsp2 = 6;
	int k = 8;
	int kp1 = 9;
	int kp2 = 10;

	std::cout << "Fortran indicies" << '\n';
	for (size_t i = nsp2; i <= kp1 ; i++)
	{
		int limit2 = kp2 - i;
		for (size_t j = 1; j <= limit2; j++)
		{
			std::cout << i << ' ' << j << '\n';
		}
	}

	std::cout << "c++ indicies" << '\n';

	for (size_t i = nsp1; i < kp1; i++)
	{
		size_t limit2 = kp1 - i;
		for (size_t j = 0; j < limit2; j++)
		{
			std::cout << i << ' ' << j << '\n';
		}
	}

	return 0;
}