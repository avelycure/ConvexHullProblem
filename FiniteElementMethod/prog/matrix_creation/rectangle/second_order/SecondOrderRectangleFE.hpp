#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "../../../common/init/InitFuncs.hpp"
#include "../../../common/classes/mesh/Point.hpp"
#include "../../border_conditions/BorderConditions.hpp"
#include "../../../common/classes/contribution_matrix/rectangle/second_order/SecondOrderRectangleRightPart.hpp"
#include "../../../common/classes/contribution_matrix/rectangle/second_order/SecondOrderRectangleContributionMatrix.hpp"
#include "../../../common/classes/system/SystemParameters.hpp"

void solveWithSecondOrderRectangleFE(SecondOrderRectangleContributionMatrix *&contributionMatrix,
                                     SecondOrderRectangleRightPart *&localRigthParts,
                                     Point **&coordinateMesh,
                                     double **&matrixPressure,
                                     double *&rightPart,
                                     SystemParameters &systemParameters);

void createLocalMatrixForEveryRectangleElementSecondOrder(SecondOrderRectangleContributionMatrix *&contributionMatrixParam,
                                                          Point **&coordinateMeshParam,
                                                          SecondOrderRectangleRightPart *&rightPartParam,
                                                          SystemParameters &systemParameters);

void createLocalContributionMatrixForRectangleElementSecondOrder(SecondOrderRectangleContributionMatrix &localMatrix,
                                                      Point pointI,
                                                      Point pointJ,
                                                      Point pointK,
                                                      Point pointL,
                                                      Point pointM,
                                                      Point pointN,
                                                      Point pointR,
                                                      Point pointQ,
                                                      SecondOrderRectangleRightPart &localRightPart,
                                                      SystemParameters &systemParameters);