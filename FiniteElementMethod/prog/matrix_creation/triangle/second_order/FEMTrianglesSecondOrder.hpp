#include <cmath>
#include <string>
#include <fstream>
#include <iostream>
#include "../../../common/init/InitFuncs.hpp"
#include "../common/TriangleCommon.hpp"
#include "../../../common/classes/mesh/Point.hpp"
#include "../../../common/classes/system/SystemParameters.hpp"
#include "../../border_conditions/BorderConditions.hpp"
#include "../../../common/classes/contribution_matrix/triangle/TriangleRightPart.hpp"
#include "../../../common/classes/contribution_matrix/triangle/TriangleContributionMatrix.hpp"
#include "../../../common/classes/contribution_matrix/triangle/TriangleRightPartSecondOrder.hpp"
#include "../../../common/classes/contribution_matrix/triangle/TriangleContributionMatrixSecondOrder.hpp"

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