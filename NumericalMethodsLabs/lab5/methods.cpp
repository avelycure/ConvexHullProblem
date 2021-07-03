#include "header.hpp"

void solveWithQuadratureMethod(double a, double b, int N)
{
	double h = (b - a) / (N - 1.0);

	double **A = new double *[N];
	for (int i = 0; i < N; i++)
	{
		A[i] = new double[N];
	}

	double *x = new double[N];
	double *d = new double[N];

	for (int i = 0; i < N; i++)
	{
		x[i] = a + i * h;
		d[i] = f(x[i]);
	}

	for (int i = 0; i < N; i++)
	{
		A[i][0] = -h / 2.0 * K(x[i], x[0]);
		for (int k = 1; k < N - 1; k++)
		{
			A[i][k] = -h * K(x[i], x[k]);
		}
		A[i][N - 1] = -h / 2.0 * K(x[i], x[N - 1]);

		A[i][i] += 1.0; //in eq there we add function twice in diagonal
	}

	findx(N, d, A);
	backwordGauss(N, d, A);

	cout << "Error in quadratic method(1 example): " << computeErrorInFirstExample(N, d) << endl;

	exportSolutionToFile(N, x, d);

	for (int i = 0; i < N; i++)
	{
		delete[] A[i];
	}
	delete[] A;
	delete[] d;
	delete[] x;
}

void solveWithMethodOfSimpleIterations(double a, double b, int N)
{
	double h = (b - a) / (N - 1.);
	double I;
	double max, temp;
	double eps = 1e-9;
	int k = 0;

	double *x = new double[N];
	double *uk = new double[N];
	double *uk1 = new double[N];
	for (int i = 0; i < N; i++)
	{
		x[i] = a + i * h;
		uk[i] = f(x[i]);
	}

	do
	{
		max = 0.0;
		temp = 0.0;

		for (int i = 0; i < N; i++)
		{
			I = 0.0;
			for (int j = 0; j < N - 1; j++)
			{
				I += 0.5 * h * (K(x[i], x[j]) * uk[j] + K(x[i], x[j + 1]) * uk[j + 1]);
			}
			uk1[i] = f(x[i]) + I;

			temp = fabs(uk1[i] - uk[i]);
			if (max < temp)
				max = temp;
		}

		copy(N, uk, uk1);
		k++;

	} while (max > eps);

	cout << k << " iterations" << endl;

	cout << "Error in method of simle iterations(1 example): " << computeErrorInFirstExample(N, uk1) << endl;
	exportSolutionToFile(N, x, uk1);

	delete[] uk1;
	delete[] uk;
	delete[] x;
}

void solveEquationWithDegenerateKernel(double a, double b, int N, int m)
{
	double h = (b - a) / (N - 1.0);

	double **alfa = new double *[m];
	for (int i = 0; i < m; i++)
		alfa[i] = new double[m];

	double *betta = new double[m]; //right part

	double *x = new double[N]; //number of nodes while computing integral
	double *u = new double[N]; //solution

	for (int i = 0; i < N; i++)
		x[i] = a + i * h;

	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < m; j++)
		{
			alfa[i][j] = 0.0;
			for (int k = 0; k < N - 1; k++)
				alfa[i][j] -= 0.5 * h * (psi(x[k], i) * phi(x[k], j) + psi(x[k + 1], i) * phi(x[k + 1], j));
		}

		alfa[i][i] += 1.0;
		betta[i] = 0.0;
		for (int k = 0; k < N - 1; k++)
			betta[i] += 0.5 * h * (psi(x[k], i) * f(x[k]) + psi(x[k + 1], i) * f(x[k + 1]));
	}

	findx(m, betta, alfa);
	backwordGauss(m, betta, alfa);

	for (int i = 0; i < N; i++)
	{
		u[i] = f(x[i]);
		for (int j = 0; j < m; j++)
			u[i] += betta[j] * phi(x[i], j);
	}

	cout << "Error while solving equation with degenerate kernel(1 example): " << computeErrorInFirstExample(N, u) << endl;

	exportSolutionToFile(N, x, u);

	for (int i = 0; i < m; i++)
		delete[] alfa[i];
	delete[] alfa;
	delete[] betta;
	delete[] x;
	delete[] u;
}

void solveEquationWithDegenerateKernelChebyshev(double a, double b, int N, int m) //m - number of functions fi and psi
{
	double h = (b - a) / (N - 1.);

	double **alfa = new double *[m];
	for (int i = 0; i < m; i++)
		alfa[i] = new double[m];

	double *betta = new double[m]; //right part

	double *x = new double[N]; //number of nodes while compting integral
	double *u = new double[N]; //solution

	for (int i = 0; i < N; i++)
		x[i] = a + i * h;

	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < m; j++)
		{
			alfa[i][j] = 0.;
			for (int k = 0; k < N - 1; k++)
				alfa[i][j] -= 0.5 * h * (psiC(x[k], i) * phiC(x[k], j) + psiC(x[k + 1], i) * phiC(x[k + 1], j));
		}

		alfa[i][i] += 1.;
		betta[i] = 0.;
		for (int k = 0; k < N - 1; k++)
			betta[i] += 0.5 * h * (psiC(x[k], i) * f(x[k]) + psiC(x[k + 1], i) * f(x[k + 1]));
	}

	findx(m, betta, alfa);
	backwordGauss(m, betta, alfa);

	for (int i = 0; i < N; i++)
	{
		u[i] = f(x[i]);
		for (int j = 0; j < m; j++)
			u[i] += betta[j] * phi(x[i], j);
	}

	cout << "Error while solving equation with degenerate kernel(Chebyshev)(1 example): " << computeErrorInFirstExample(N, u) << endl;

	exportSolutionToFile(N, x, u);

	for (int i = 0; i < m; i++)
		delete[] alfa[i];
	delete[] alfa;
	delete[] betta;
	delete[] x;
	delete[] u;
}

void singularKernel(int N)
{
	double **A = new double *[N + 1];
	for (int i = 0; i < N + 1; i++)
	{
		A[i] = new double[N + 1];
	}
	double *d = new double[N + 1]; //right part

	double h = 2.0 * PI / N;

	double *cx = new double[N];
	double *cy = new double[N];
	double *kx = new double[N]; //=nx
	double *ky = new double[N]; //=ny

	for (int i = 0; i < N; i++)
	{
		cx[i] = cos(h * i);
		cy[i] = sin(h * i);
		kx[i] = cos(h * (i + 0.5));
		ky[i] = sin(h * (i + 0.5));
		d[i] = f(h * (i + 0.5));
	}

	d[N] = 0.0;
	A[N][N] = 0.0;

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{ //nxQx+nyQy
			A[i][j] = 1.0 / N * (ky[i] * (kx[i] - cx[j]) - kx[i] * (ky[i] - cy[j])) / ((kx[i] - cx[j]) * (kx[i] - cx[j]) + (ky[i] - cy[j]) * (ky[i] - cy[j]));
		}
		A[i][N] = 1.0; //extra column
		A[N][i] = 1.0; //extra row
	}

	findx(N + 1, d, A);
	backwordGauss(N + 1, d, A);

	double temp = 0.0;
	double err = 0.0;

	for (int i = 0; i < N; i++)
	{
		temp = fabs(d[i] - 2 * sin(7 * i * h));
		if (err < temp)
			err = temp;
	}

	cout << "Error while solving equation with singular kernel: " << err << endl;

	exportSolutionToFile(N, cx, cy, d);

	for (int i = 0; i < N + 1; i++)
	{
		delete[] A[i];
	}

	delete[] A;
	delete[] d;
	delete[] cx;
	delete[] cy;
	delete[] kx;
	delete[] ky;
}

double computeErrorInFirstExample(int N, double *y)
{
	double temp = 0.0;
	double err = 0.0;

	for (int i = 0; i < N; i++)
	{
		temp = fabs(y[i] - 1.0);
		if (err < temp)
		{
			err = temp;
		}
	}

	return err;
}

bool findx(const int DIM, double *b, double **a) //straight gauss with checking degeneracy
{
	double max;
	double *temp1;
	double temp2;
	bool flag = true;

	for (int k = 0, h; k < DIM; k++)
	{
		max = fabs(a[k][k]); //writing a diagonal element to a comparison variable
		temp1 = a[k];		 //writing a reference to the string of a diagonal element (if the condition in the loop is never met)
		h = k;				 //writing the line number (for the reason above)

		if (fabs(a[k][k]) < EPSILON)
		{
			flag = true;
		} //if the first element is zero, then activate the flag
		else
		{
			flag = false;
		} //if not zero, then deactivate the flag

		for (int i = (k + 1); i < DIM; i++)
		{
			if (fabs(a[i][k]) > max && fabs(a[i][k]) > EPSILON) //if the element is greater than the previous max in the column and not zero
			{
				max = fabs(a[i][k]);
				temp1 = a[i];
				h = i;
				flag = false;
			}; //remember the line, change the max and deactivate the flag
		}

		if (!flag) //if the column is not null
		{
			a[h] = a[k]; //string swap
			a[k] = temp1;
			temp2 = b[h]; //swapping column elements
			b[h] = b[k];
			b[k] = temp2;
		}
		else
		{
			cout << "Matrix is degenerate\n";
			return false;
		}

		for (int i = k + 1; i < DIM; i++)
		{
			b[i] = b[i] - b[k] / a[k][k] * a[i][k]; //substraction from column elements
			for (int j = (DIM - 1); j >= 0; j--)
			{
				a[i][j] = a[i][j] - a[k][j] / a[k][k] * a[i][k]; //vanishing of all elements under the diagonal
			}
		}
	}
	return true;
}

void backwordGauss(const int DIM, double *b, double **a) //backwords gauss
{
	for (int i = (DIM - 1); i >= 0; i--)
	{
		for (int j = i + 1; j < DIM; j++)
		{
			b[i] = b[i] - a[i][j] * b[j]; //subtraction from the column of the right side of all elements of the matrix of the same row, except for the element on the global diagonal
		}
		b[i] = b[i] / a[i][i]; //dividing the column element on the right side by the main diagonal element
	}
}

double singularKernelNewVersion(int N)
{
	double **A = new double *[N + 1];
	for (int i = 0; i < N + 1; i++)
	{
		A[i] = new double[N + 1];
	}
	double *d = new double[N + 1]; //right part

	double h = 2.0 * PI / N;

	double *cx = new double[N + 1];
	double *cy = new double[N + 1];
	double *kx = new double[N + 1]; //=nx
	double *ky = new double[N + 1]; //=ny

	for (int i = 0; i < N + 1; i++)
	{
		cx[i] = cos(h * i);
		cy[i] = sin(h * i);
		kx[i] = cos(h * (i + 0.5));
		ky[i] = sin(h * (i + 0.5));
		d[i] = f(h * (i + 0.5));
	}

	//d[N] = 0.0;
	//A[N][N] = 0.0;

	for (int i = 0; i < N + 1; i++)
	{
		for (int j = 0; j < N + 1; j++)
		{ //nxQx+nyQy
			A[i][j] = 1.0 / N * (ky[i] * (kx[i] - cx[j]) - kx[i] * (ky[i] - cy[j])) / ((kx[i] - cx[j]) * (kx[i] - cx[j]) + (ky[i] - cy[j]) * (ky[i] - cy[j]));
		}
	}

	findx(N + 1, d, A);
	backwordGauss(N + 1, d, A);

	double temp = 0.0;
	double err = 0.0;

	for (int i = 0; i < N + 1; i++)
	{
		temp = fabs(d[i] - 2 * sin(VARIANT_NUMBER * i * h));
		if (err < temp)
			err = temp;
	}

	cout << "Error while solving equation with singular kernel: " << err << endl;
	double R = d[N];
	cout << "Value of R: " << R << endl;

	exportSolutionToFile(N + 1, cx, cy, d);

	for (int i = 0; i < N + 1; i++)
	{
		delete[] A[i];
	}

	delete[] A;
	delete[] d;
	delete[] cx;
	delete[] cy;
	delete[] kx;
	delete[] ky;
	return R;
}