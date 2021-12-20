#pragma once
#include <cmath>
#include <string>
#include <fstream>
#include <iostream>
#include "../../common/init/InitFuncs.hpp"
#include "../../common/classes/mesh/Point.hpp"
#include "../../common/classes/system/SystemParameters.hpp"
#include "../border_conditions/BorderConditions.hpp"
#include "../../common/classes/contribution_matrix/triangle/TriangleRightPart.hpp"
#include "../../common/classes/contribution_matrix/triangle/TriangleContributionMatrix.hpp"
#include "../../common/classes/contribution_matrix/triangle/TriangleRightPartSecondOrder.hpp"
#include "../../common/classes/contribution_matrix/triangle/TriangleContributionMatrixSecondOrder.hpp"

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

void createGlobalPressureMatrixHConst(double **&matrixPressure,
                                      TriangleContributionMatrix *&contributionMatrix,
                                      int n);

void solveWithHConst(TriangleContributionMatrix *&contributionMatrix,
                     Point **&coordinateMesh,
                     double **&matrixPressure,
                     double *&rightPart,
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

//Quadratic
void solveWithTrianglesSecondOrder(TriangleContributionMatrixSecondOrder *&contributionMatrix,
                                   TriangleRightPartSecondOrder *&localRigthParts,
                                   Point **&coordinateMesh,
                                   double **&matrixPressure,
                                   double *&rightPart,
                                   SystemParameters &systemParameters);

void createLocalContributionMatrixForQuardaticTriangle(TriangleContributionMatrixSecondOrder localMatrix,
                                                          Point pointI,
                                                          Point pointJ,
                                                          Point pointK,
                                                          Point pointL,
                                                          Point pointM,
                                                          Point pointN,
                                                          TriangleRightPartSecondOrder localRightPart,
                                                          SystemParameters &systemParameters);

void addBorderConditionsQuadraticTriangles(double **&matrixResult,
                                           double *&rightPartParam,
                                           int n,
                                           int MATRIX_PRESSURE_SIZE,
                                           int OTHER_BORDER,
                                           int DOWN_BORDER);

void createGlobalPressureMatrixQuadraticTriangles(double **&matrixPressure,
                                                  TriangleContributionMatrixSecondOrder *&contributionMatrix,
                                                  double *&rightPartParam,
                                                  TriangleRightPartSecondOrder *&localRightPartsParam,
                                                  int n);

void createLocalMatrixForEveryElementQuadraticTriangles(TriangleContributionMatrixSecondOrder *&contributionMatrixParam,
                                                        Point **&coordinateMeshParam,
                                                        TriangleRightPartSecondOrder *&rightPartParam,
                                                        SystemParameters &systemParameters);