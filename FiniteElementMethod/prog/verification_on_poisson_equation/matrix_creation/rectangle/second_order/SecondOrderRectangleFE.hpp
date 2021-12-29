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

void createGlobalPressureMatrixForRectangleElementSecondOrder(double **&matrixPressure,
                                                              SecondOrderRectangleContributionMatrix *&contributionMatrix,
                                                              double *&rightPartParam,
                                                              SecondOrderRectangleRightPart *&localRightPartsParam,
                                                              int n);

void setCoefficients(double &c1, double &c2, double &c3, double &c4,
                     double &c5, double &c6, double &c7, double &c8,
                     double &c9, double &c10, double &c11, double &c12,
                     double &c13, double &c14, double &c15, double &c16,
                     double &c17, double &c18, double &c19, double &c20,
                     double &c21, double &c22, double &c23, double &c24,
                     double &c25, double &c26, double &c27, double &c28,
                     double &c29, double &c30, double &c31, double &c32,
                     double &c33, double &c34, double &c35, double &c36,
                     double &c37, double &c38, double &c39, double &c40, double &c41,
                     double A1, double A2, double A3, double A4,
                     double s1, double s2, double s3, double s4, double s5, double s6, double s7, double s8,
                     double z1, double z2, double z3, double z4, double z5);

void setFormFunctionsCoefficients(double *&a,
                                  double *&b,
                                  double *&c,
                                  double *&d,
                                  double *&e,
                                  double *&f,
                                  double *&g,
                                  double *&t,
                                  Point pointI,
                                  Point pointM);