#include "header.hpp"
/*
	на 3, 7 примерах использован шаг 0.01; на остальных 0.1
*/
int main()
{
	//integration limits
	double a = 0.1;
	double b = 1.0;
	int teylorCoeffsNum = 5;

	double h = 0.01;
	int N = round((b - a) / h) + 1;

	//solveWithQuadratureMethod(a, b, N);
	//solveWithMethodOfSimpleIterations(a, b, N);
	//solveEquationWithDegenerateKernel(a, b, N, teylorCoeffsNum);

	//part 2
	//N = 500;
	//singularKernelNewVersion(N);
	//singularKernel(N);

	plotDependencyRfromNInSingularEquation();

	return 0;
}

void plotDependencyRfromNInSingularEquation()
{
	ofstream fout;
	fout.open(FILE_DEPENDENCY_NAME, ios_base::out | ios_base::trunc);
	int pointsNum = 1000;

	for (int i = 50; i < pointsNum; i += 50)
	{
		fout << i << " " << singularKernelNewVersion(i) << endl;
	}
	fout.close();
	fout.clear();
}
