#include <string>
#include <fstream>
#include <iostream>
#include <cmath>
#include "../common/classes/system/SystemParameters.hpp"

double compareWithAnalyticNormSum(double *&solutionFE,
                                  int n,
                                  SystemParameters systemParameters);

double compareWithAnalyticNormMax(double *&solutionFE,
                                  int n,
                                  SystemParameters systemParameters);

double analytic(double x,
                double y);

bool readSolution(const std::string fileNameVector,
                  double *&solution,
                  const int &n);