#pragma once
#include <fstream>
#include<iostream>
#include<iomanip>
#include<limits>
#include<cmath>

void addBorderConditionsOnPressureValues(double **&matrixResult,
                                         double *&rightPart,
                                         int n,
                                         int MATRIX_PRESSURE_SIZE,
                                         double LOW_BORDER,
                                         double HIGH_BORDER);

void addBorderConditionsPressureDerivatives(double **&matrixResult,
                                            int n,
                                            double h,
                                            int MATRIX_PRESSURE_SIZE,
                                            double LOW_BORDER,
                                            double HIGH_BORDER);

void addBorderConditionsForRectangleSecondOrder(double **&matrixResult,
                                                double *&rightPart,
                                                int n,
                                                int MATRIX_PRESSURE_SIZE,
                                                double LOW_BORDER,
                                                double HIGH_BORDER);