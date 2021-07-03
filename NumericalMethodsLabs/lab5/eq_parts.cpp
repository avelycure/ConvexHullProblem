#include "header.hpp"

//right part
double f(double x)
{
	//part1
	//return 0.5 * (1.0 + sin(x));
	//return x * x + sqrt(x);

	//part2
	return cos((VARIANT_NUMBER / 2.0) * x);
	//return sin(((VARIANT_NUMBER + 1) / 2.0) * x);
}

//kernel
double K(double x, double s)
{
	return 0.5 * (1.0 - x * cos(x * s));
}
/*
	funcs for solving equations with degenerate kernel
	phi and psi are systems of independent functions
*/
double phi(double x, int i)
{
	if (i == 0)
		return 0.5 * (1.0 - x);

	if (i == 1)
		return 0.25 * pow(x, 3);

	if (i == 2)
		return -1.0 / 48.0 * pow(x, 5);

	if (i == 3)
		return 1.0 / 1440.0 * pow(x, 7);

	if (i == 4)
		return (-1.0 / 80640.0) * pow(x, 9);
}

double psi(double s, int i)
{
	if (i == 0)
		return 1.0;

	if (i == 1)
		return pow(s, 2);

	if (i == 2)
		return pow(s, 4);

	if (i == 3)
		return pow(s, 6);

	if (i == 4)
		return pow(s, 8);
}

double phiC(double x, int i)
{
	if (i == 0)
		return 0.5 * (1.0 - x);

	if (i == 1)
		return 0.25 * pow(x, 3);

	if (i == 2)
		return -0.0208332 * pow(x, 5);

	if (i == 3)
		return 0.000694148 * pow(x, 7);
}

double psiC(double s, int i)
{
	if (i == 0)
		return 1.0;

	if (i == 1)
		return pow(s, 2);

	if (i == 2)
		return pow(s, 4);

	if (i == 3)
		return pow(s, 6);
}
