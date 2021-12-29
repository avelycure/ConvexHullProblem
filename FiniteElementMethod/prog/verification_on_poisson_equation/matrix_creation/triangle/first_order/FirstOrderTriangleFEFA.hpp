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

//Functional approach functions
void solveWithFirstOrderTriangleFEConstantHeight(TriangleContributionMatrix *&contributionMatrix,
                                                 Point **&coordinateMesh,
                                                 double **&matrixPressure,
                                                 double *&rightPart,
                                                 SystemParameters &systemParameters);

void createLocalContributionMatrixConstantHeight(TriangleContributionMatrix localMatrix,
                                                    Point pointI,
                                                    Point pointJ,
                                                    Point pointK);

void createLocalMatrixForEveryElementConstantHeight(TriangleContributionMatrix *&contributionMatrixParam,
                                                    Point **&coordinateMeshParam,
                                                    int n);

void createGlobalPressureMatrixConstantHeight(double **&matrixPressure,
                                              TriangleContributionMatrix *&contributionMatrix,
                                              int n);