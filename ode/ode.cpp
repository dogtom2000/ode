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

	std::cout << out;
	return 0;
}