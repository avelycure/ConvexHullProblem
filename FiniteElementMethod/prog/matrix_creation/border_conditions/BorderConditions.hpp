#pragma once
#include <fstream>

void addBorderConditionsOnPressureValues(double **&matrixResult,
                                         double *&rightPart,
                                         int n,
                                         int MATRIX_PRESSURE_SIZE,
                                         double LOW_BORDER,
                                         double HIGH_BORDER);