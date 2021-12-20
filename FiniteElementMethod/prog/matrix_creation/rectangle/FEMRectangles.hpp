#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "../../common/init/InitFuncs.hpp"
#include "../../common/classes/mesh/Point.hpp"
#include "../border_conditions/BorderConditions.hpp"
#include "../../common/classes/contribution_matrix/rectangle/RectangleContributionMatrix.hpp"
#include "../../common/classes/contribution_matrix/rectangle/RectangleRightPart.hpp"
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

void createGlobalPressureMatrixForRectangleElement(double **&matrixPressure,
                                                   RectangleContributionMatrix *&contributionMatrix,
                                                   double *&rightPartParam,
                                                   RectangleRightPart *&localRightPartsParam,
                                                   int n);