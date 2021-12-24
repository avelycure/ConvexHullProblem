#include <string>
#include <fstream>
#include <iostream>
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

double analyticOfLaplasEquation(double x, double y, SystemParameters systemParameters);

bool readSolution(const std::string fileNameVector,
                  double *&solution,
                  const int &n);

void outputLaplasSolution(Point **&coordinateMesh,
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

double a(int n);

double b(int n);

double c(int n);

double d(int n);