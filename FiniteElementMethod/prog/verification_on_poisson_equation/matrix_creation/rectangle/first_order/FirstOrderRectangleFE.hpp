#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "../../../common/init/InitFuncs.hpp"
#include "../../../common/classes/mesh/Point.hpp"
#include "../../border_conditions/BorderConditions.hpp"
#include "../../../common/classes/contribution_matrix/rectangle/first_order/FirstOrderRectangleContributionMatrix.hpp"
#include "../../../common/classes/contribution_matrix/rectangle/first_order/FirstOrderRectangleRightPart.hpp"
#include "../../../common/classes/system/SystemParameters.hpp"

//Rectangle
void createLocalMatrixForEveryRectangleElement(FirstOrderRectangleContributionMatrix *&contributionMatrixParam,
                                               Point **&coordinateMeshParam,
                                               FirstOrderRectangleRightPart *&rightPartParam,
                                               SystemParameters &systemParameters);

void createLocalContributionMatrixForRectangleElement(FirstOrderRectangleContributionMatrix &localMatrix,
                                                      Point pointI,
                                                      Point pointJ,
                                                      Point pointK,
                                                      Point pointM,
                                                      FirstOrderRectangleRightPart &localRightPart,
                                                      SystemParameters &systemParameters);

void solveWithFirstOrderRectangleFE(FirstOrderRectangleContributionMatrix *&contributionMatrix,
                                      FirstOrderRectangleRightPart *&localRigthParts,
                                      Point **&coordinateMesh,
                                      double **&matrixPressure,
                                      double *&rightPart,
                                      SystemParameters &systemParameters);

void createGlobalPressureMatrixForRectangleElement(double **&matrixPressure,
                                                   FirstOrderRectangleContributionMatrix *&contributionMatrix,
                                                   double *&rightPartParam,
                                                   FirstOrderRectangleRightPart *&localRightPartsParam,
                                                   int n);