#include "header.hpp"

void copy(const int DIM, double **a, double **b)
{
	for (int i = 0; i < DIM; i++)
	{
		for (int j = 0; j < DIM; j++)
		{
			a[i][j] = b[i][j];
		}
	}
}

void copy(const int DIM, double *a, double *b)
{
	for (int j = 0; j < DIM; j++)
	{
		a[j] = b[j];
	}
}

void exportSolutionToFile(int N, double *x, double *Y)
{
	ofstream fout;
	fout.open(FILE_RESULT_NAME, ios_base::out | ios_base::trunc);

	for (int i = 0; i < N; i++)
	{
		fout << x[i] << " " << Y[i] << endl;
	}
	fout.close();
	fout.clear();
};

void exportSolutionToFile(int N, double *x, double *y, double *g)
{
	ofstream fout;
	fout.open(FILE_RESULT_SINGLE_NAME, ios_base::out | ios_base::trunc);

	for (int i = 0; i < N; i++)
	{
		fout << x[i] << " " << y[i] << " " << g[i] << endl;
	}
	fout.close();
	fout.clear();
};

void displayMatrix(const int DIM, double **const a)
{
	std::cout << endl;
	for (int i = 0; i < DIM; i++)
	{
		for (int j = 0; j < DIM; j++)
		{
			std::cout << a[i][j] << " ";
		}
		std::cout << endl;
	}
	std::cout << endl;
}

void displayVector(const int DIM, const double *const b)
{
	std::cout << '(';
	for (int i = 0; i < DIM; i++)
	{
		cout << b[i] << ",";
	}
	std::cout << ')';
}
