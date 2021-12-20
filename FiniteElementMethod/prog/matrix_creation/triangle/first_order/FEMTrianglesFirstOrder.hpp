#pragma once
#include <cmath>
#include <string>
#include <fstream>
#include <iostream>
#include "../common/TriangleCommon.hpp"
#include "../../../common/init/InitFuncs.hpp"
#include "../../../common/classes/mesh/Point.hpp"
#include "../../../common/classes/system/SystemParameters.hpp"
#include "../../border_conditions/BorderConditions.hpp"
#include "../../../common/classes/contribution_matrix/triangle/TriangleRightPart.hpp"
#include "../../../common/classes/contribution_matrix/triangle/TriangleContributionMatrix.hpp"
#include "../../../common/classes/contribution_matrix/triangle/TriangleRightPartSecondOrder.hpp"
#include "../../../common/classes/contribution_matrix/triangle/TriangleContributionMatrixSecondOrder.hpp"

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

void solveWithHLinear(TriangleContributionMatrix *&contributionMatrix,
                     TriangleRightPart *&localRigthParts,
                     Point **&coordinateMesh,
                     double **&matrixPressure,
                     double *&rightPart,
                     SystemParameters &systemParameters);