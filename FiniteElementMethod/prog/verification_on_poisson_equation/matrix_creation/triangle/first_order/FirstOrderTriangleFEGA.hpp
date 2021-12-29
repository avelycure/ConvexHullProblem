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

//Galerkin approach functions

void solveWithFirstOrderTriangleFE(TriangleContributionMatrix *&contributionMatrix,
                                   TriangleRightPart *&localRigthParts,
                                   Point **&coordinateMesh,
                                   double **&matrixPressure,
                                   double *&rightPart,
                                   SystemParameters &systemParameters);

void createLocalContributionMatrix(TriangleContributionMatrix localMatrix,
                                   Point pointI,
                                   Point pointJ,
                                   Point pointK,
                                   TriangleRightPart localRightPart,
                                   SystemParameters &systemParameters);

void createLocalMatrixes(TriangleContributionMatrix *&contributionMatrixParam,
                         Point **&coordinateMeshParam,
                         TriangleRightPart *&rightPartParam,
                         SystemParameters &systemParameters);

void createGlobalPressureMatrix(double **&matrixPressure,
                                TriangleContributionMatrix *&contributionMatrix,
                                double *&rightPartParam,
                                TriangleRightPart *&localRightPartsParam,
                                int n);