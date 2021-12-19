#pragma once
#include <cmath>
#include <string>
#include <fstream>
#include <iostream>
#include "../../common/init/InitFuncs.hpp"
#include "../../common/classes/mesh/Point.hpp"
#include "../../common/classes/system/SystemParameters.hpp"
#include "../../common/classes/contribution_matrix/triangle/TriangleRightPart.hpp"
#include "../../common/classes/contribution_matrix/triangle/TriangleContributionMatrix.hpp"

//Finite elements method funcs
double countArea(Point pointI,
                 Point pointJ,
                 Point pointK);

void createLocalContributionMatrixForHConst(TriangleContributionMatrix localMatrix,
                                            Point pointI,
                                            Point pointJ,
                                            Point pointK);

void createLocalMatrixForEveryElementHConst(TriangleContributionMatrix *&contributionMatrixParam,
                                            Point **&coordinateMeshParam,
                                            int n);

void addBorderConditionsHConst(double **&matrixResult,
                               int n,
                               int MATRIX_PRESSURE_SIZE,
                               int OTHER_BORDER,
                               int DOWN_BORDER);

void createGlobalPressureMatrixHConst(double **&matrixPressure,
                                      TriangleContributionMatrix *&contributionMatrix,
                                      int n);

void solveWithHConst(TriangleContributionMatrix *&contributionMatrix,
                     Point **&coordinateMesh,
                     double **&matrixPressure,
                     SystemParameters &systemParameters);

int createLocalContributionMatrixForHLinearTop(TriangleContributionMatrix localMatrix,
                                               Point pointI,
                                               Point pointJ,
                                               Point pointK,
                                               TriangleRightPart localRightPart,
                                               SystemParameters &systemParameters);

int createLocalContributionMatrixForHLinearBottom(TriangleContributionMatrix localMatrix,
                                                  Point pointI,
                                                  Point pointJ,
                                                  Point pointK,
                                                  TriangleRightPart localRightPart,
                                                  SystemParameters &systemParameters);
void createLocalMatrixForEveryElementHLinear(TriangleContributionMatrix *&contributionMatrixParam,
                                             Point **&coordinateMeshParam,
                                             TriangleRightPart *&rightPartParam,
                                             SystemParameters &systemParameters);

void createGlobalPressureMatrixHLinear(double **&matrixPressure,
                                       TriangleContributionMatrix *&contributionMatrix,
                                       double *&rightPartParam,
                                       TriangleRightPart *&localRightPartsParam,
                                       int n);

void addBorderConditionsHLinear(double **&matrixResult,
                                double *&rightPartParam, int n,
                                int MATRIX_PRESSURE_SIZE,
                                int OTHER_BORDER,
                                int DOWN_BORDER);

/**
 * Right border == left border
 * */
void addBorderConditionsToLeftAndRight(double **&matrixResult,
                                       int n,
                                       double h,
                                       int MATRIX_PRESSURE_SIZE,
                                       double TOP_BORDER,
                                       double BOTTOM_BORDER);

int solveWithHLinearWithDerBC(TriangleContributionMatrix *&contributionMatrix,
                              TriangleRightPart *&localRigthParts,
                              Point **&coordinateMesh,
                              double **&matrixPressure,
                              double *&rightPart,
                              SystemParameters &systemParameters);

int solveWithHConstBCLR(TriangleContributionMatrix *&contributionMatrix,
                        Point **&coordinateMesh,
                        double **&matrixPressure,
                        SystemParameters &systemParameters);

int solveWithHLinear(TriangleContributionMatrix *&contributionMatrix,
                     TriangleRightPart *&localRigthParts,
                     Point **&coordinateMesh,
                     double **&matrixPressure,
                     double *&rightPart,
                     SystemParameters &systemParameters);