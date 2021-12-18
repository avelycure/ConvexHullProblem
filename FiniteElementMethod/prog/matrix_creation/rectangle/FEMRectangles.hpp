#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "../../common/init/InitFuncs.hpp"
#include "../../common/classes/mesh/Point.hpp"
#include "../../common/classes/contribution_matrix/RectangleContributionMatrix.hpp"
#include "../../common/classes/contribution_matrix/RectangleRightPart.hpp"
#include "../../common/classes/system/SystemParameters.hpp"

//Rectangle
void createLocalMatrixForEveryRectangleElement(RectangleContributionMatrix *&contributionMatrixParam,
                                               Point **&coordinateMeshParam,
                                               RectangleRightPart *&rightPartParam,
                                               SystemParameters &systemParameters);

void createLocalContributionMatrixForRectangleElement(RectangleContributionMatrix &localMatrix,
                                                      Point pointI,
                                                      Point pointJ,
                                                      Point pointK,
                                                      Point pointM,
                                                      RectangleRightPart &localRightPart,
                                                      SystemParameters &systemParameters);

void solveWithRectangleFiniteElements(RectangleContributionMatrix *&contributionMatrix,
                                      RectangleRightPart *&localRigthParts,
                                      Point **&coordinateMesh,
                                      double **&matrixPressure,
                                      double *&rightPart,
                                      SystemParameters &systemParameters);

void createGlobalPressureMatrixForRectangleElement(
    double **&matrixPressure, RectangleContributionMatrix *&contributionMatrix,
    double *&rightPartParam, RectangleRightPart *&localRightPartsParam, int n);

void addBorderConditionsForRectnangleElements(double **&matrixResult,
                                              double *&rightPartParam,
                                              int n,
                                              int MATRIX_PRESSURE_SIZE,
                                              int OTHER_BORDER,
                                              int DOWN_BORDER);

void addBorderConditionsForRectnangleElementsToLeftAndRight(double **&matrixResult,
                                       int n,
                                       double h,
                                       int MATRIX_PRESSURE_SIZE,
                                       double TOP_BORDER,
                                       double BOTTOM_BORDER);