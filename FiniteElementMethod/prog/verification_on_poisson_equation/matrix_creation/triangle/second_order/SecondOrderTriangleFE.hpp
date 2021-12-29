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


void solveWithSecondOrderTriangleFE(TriangleContributionMatrixSecondOrder *&contributionMatrix,
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

void createGlobalPressureMatrixQuadraticTriangles(double **&matrixPressure,
                                                  TriangleContributionMatrixSecondOrder *&contributionMatrix,
                                                  double *&rightPartParam,
                                                  TriangleRightPartSecondOrder *&localRightPartsParam,
                                                  int n);

void createLocalMatrixForEveryElementQuadraticTriangles(TriangleContributionMatrixSecondOrder *&contributionMatrixParam,
                                                        Point **&coordinateMeshParam,
                                                        TriangleRightPartSecondOrder *&rightPartParam,
                                                        SystemParameters &systemParameters);

void setCoefficients(double &c1, double &c2, double &c3, double &c4, double &c5, double &c6, double &c7,
                     double &c8, double &c9, double &c10, double &c11, double &c12, double &c13,
                     double &c14, double &c15, double &c16, double &c18, double &c19,
                     double s1, double s2, double s3, double s4, double s5, double s6, double s7,
                     double k1, double k2,
                     double A1, double A2, double A3, double A4,
                     double zi);

void setFormFunctionsCoefficients(double *&a,
                                  double *&b,
                                  double *&c,
                                  double *&d,
                                  double *&e,
                                  double *&f,
                                  Point pointI,
                                  Point pointJ,
                                  Point pointK,
                                  Point pointM,
                                  Point pointN);