#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>

using namespace std;
void solveEquation(const int size);

void allocateMemory(double **&A, double *&B, double *&X, const int &n);
void freeMemory(double **&A, double *&B, double *&X, const int &n);

bool readData(const string fileNameMatrix, const string fileNameVector, double **&matrixA, double *&vectorB, const int &n);

int WriteVector(string fileNameOutput, double *&vector, const int &n);

bool solveWithGauss(double **&A, double *&B, double *&X, const int &size);

bool matrixIsPrepared(double **&A, double *&B, const int &i, const int &size);

void diagonalizeEquation(double **&A, double *&B, double *&X, const int &size);

int searchMaxInColumn(double **&A, const int &currentMinor, const int &size);

void swapRows(double **&A, double *&B, const int &i1, const int &i2);

bool isDegenerate(double **&A, const int &i, const int &size);
