#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

const int VARIANT_NUMBER = 4;
const string FILE_RESULT_NAME = "solution.txt";
const string FILE_RESULT_SINGLE_NAME = "solutionSingle.txt";
const string FILE_DEPENDENCY_NAME = "dependency.txt";
const double EPSILON = 1e-16; //const for comparing with zero
const double PI = 3.141592;

//Equation parts
double phi(double x, int i);
double psi(double x, int i);
double phiC(double x, int i);
double psiC(double s, int i);
double phi(double x, int i);
double psi(double s, int i);
double f(double x);
double K(double x, double s);

//part 1
void solveWithQuadratureMethod(double a, double b, int N);
void solveWithMethodOfSimpleIterations(double a, double b, int N);
void solveEquationWithDegenerateKernel(double a, double b, int N, int m);
void solveEquationWithDegenerateKernelChebyshev(double a, double b, int N, int m);

//part 2
void singularKernel(int N);
double singularKernelNewVersion(int N);
void plotDependencyRfromNInSingularEquation();

void backwordGauss(const int DIM, double *b, double **a); //backword gouss
bool findx(const int DIM, double *b, double **a);         //straight gauss with out degeneracy matrix

void copy(const int DIM, double **a, double **b); //copy matrixes
void copy(const int DIM, double *a, double *b);   //copy vectors
void exportSolutionToFile(int N, double *x, double *Y);
void exportSolutionToFile(int N, double *x, double *y, double *g);

void displayMatrix(const int DIM, double **const a);
void displayVector(const int DIM, const double *const b);
double computeErrorInFirstExample(int N, double *y);
double phiC(double x, int i);
double psiC(double s, int i);