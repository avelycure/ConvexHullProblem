#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>
#include "../common/classes/system/SystemParameters.hpp"
#include "../common/classes/mesh/Point.hpp"

double compareWithAnalyticNormSum(Point **&coordinateMesh,
                                  int n,
                                  SystemParameters systemParameters);

double compareWithAnalyticNormMax(Point **&coordinateMesh,
                                  int n,
                                  SystemParameters systemParameters);

double compareWithAnalyticRectangleSecondOrderNormMax(Point **&coordinateMesh,
                                                      int n,
                                                      SystemParameters systemParameters);

double compareWithAnalyticRectangleSecondOrderNormSum(Point **&coordinateMesh,
                                                      int n,
                                                      SystemParameters systemParameters);

double analyticOfPoissonEquation(double x, double y, SystemParameters systemParameters);

bool readSolution(const std::string fileNameVector,
                  double *&solution,
                  const int &n);

void outputPoissonSolution(Point **&coordinateMesh,
                           std::string fileNameOutput,
                           const int &n,
                           SystemParameters systemParameters);

void outputDifference(double *&solution,
                      Point **&coordinateMesh,
                      std::string fileNameOutput,
                      const int &n,
                      SystemParameters systemParameters);

void outputAnalyticSolutionOnSquareMesh(Point **&coordinateMesh,
                                        std::string fileNameOutput,
                                        const int &n,
                                        SystemParameters systemParameters);

void outputDifferenceOnSquareMesh(double *&solution,
                                  Point **&coordinateMesh,
                                  std::string fileNameOutput,
                                  const int &n,
                                  SystemParameters systemParameters);