#pragma once
#include <fstream>
#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>
#include "../../common/classes/system/SystemParameters.hpp"
#include "../../common/classes/mesh/Point.hpp"

void addBorderConditionsOnPressureValues(SystemParameters systemParameters,
                                         double **&matrixResult,
                                         double *&rightPart,
                                         int n,
                                         int MATRIX_PRESSURE_SIZE);

void addBorderConditionsForRectangleSecondOrder(double **&matrixResult,
                                                double *&rightPart,
                                                int n,
                                                int MATRIX_PRESSURE_SIZE,
                                                SystemParameters systemParameters);