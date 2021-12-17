#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "common/classes/mesh/Point.hpp"
#include "common/classes/contribution_matrix/RectangleContributionMatrix.hpp"
#include "common/classes/contribution_matrix/TriangleContributionMatrix.hpp"
#include "common/classes/contribution_matrix/RectangleRightPart.hpp"
#include "common/classes/contribution_matrix/TriangleRightPart.hpp"
#include "common/classes/system/SystemParameters.hpp"
#include "libs/single_include/nlohmann/json.hpp"
#include "solvers/gauss/GaussSystemSolver.hpp"

using namespace std;

const string FILE_PARAMETERS_NAME = "data/fem_input/initial_conditions/systemParameters.json";
const string FILE_SYSTEM_NAME = "data/fem_input/initial_conditions/systemNum.json";
const string H_CONST = "hConst";
const string H_LINEAR = "hLinear";

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

void createGlobalPressureMatrixForRectangleElement(
    double **&matrixPressure, RectangleContributionMatrix *&contributionMatrix,
    double *&rightPartParam, RectangleRightPart *&localRightPartsParam, int n);

void addBorderConditionsForRectnangleElements(double **&matrixResult,
                                              double *&rightPartParam,
                                              int n,
                                              int MATRIX_PRESSURE_SIZE,
                                              int OTHER_BORDER,
                                              int DOWN_BORDER);

void addBorderConditionsForRectnangleElementsToLeftAndRight(double **&matrixResult,
                                       int n,
                                       double h,
                                       int MATRIX_PRESSURE_SIZE,
                                       double TOP_BORDER,
                                       double BOTTOM_BORDER);

//Basic funcs
void initMatrix(double **&matrix, int row, int column);
void displayMatrix(double **matrix, int row, int column);
void displayMesh(Point **coordinateMesh, int n);
void displayAllLocalMatrixes(TriangleContributionMatrix *&ContributionMatrixParam, int n);
void outputPressureMatrix(double **matrixPressure, int MATRIX_PRESSURE_SIZE);
void displayVector(double *mVector, int n);
void initVector(double *&p, int n);
void initMesh(Point **&coordinateMesh, SystemParameters &systemParameters);
void initRightPart(TriangleRightPart *&localRigthParts, int MATRIX_CONTRIBUTION_SIZE);
void initContributionMatrix(TriangleContributionMatrix *&contributionMatrix, int MATRIX_CONTRIBUTION_SIZE);
void readSystemParameters(SystemParameters &systemParameters, string &method);

//Finite elements method funcs
double countArea(Point pointI, Point pointJ, Point pointK);
void createLocalContributionMatrixForHConst(TriangleContributionMatrix localMatrix, Point pointI, Point pointJ, Point pointK);
void createLocalMatrixForEveryElementHConst(TriangleContributionMatrix *&contributionMatrixParam, Point **&coordinateMeshParam, int n);
void addBorderConditionsHConst(double **&matrixResult, int n, int MATRIX_PRESSURE_SIZE, int OTHER_BORDER, int DOWN_BORDER);
void createGlobalPressureMatrixHConst(double **&matrixPressure, TriangleContributionMatrix *&contributionMatrix, int n);
int solveWithHConst(TriangleContributionMatrix *&contributionMatrix,
                    Point **&coordinateMesh,
                    double **&matrixPressure,
                    SystemParameters &systemParameters);

int createLocalContributionMatrixForHLinearTop(TriangleContributionMatrix localMatrix,
                                               Point pointI, Point pointJ, Point pointK,
                                               TriangleRightPart localRightPart, SystemParameters &systemParameters);
int createLocalContributionMatrixForHLinearBottom(TriangleContributionMatrix localMatrix,
                                                  Point pointI, Point pointJ, Point pointK,
                                                  TriangleRightPart localRightPart, SystemParameters &systemParameters);
void createLocalMatrixForEveryElementHLinear(TriangleContributionMatrix *&contributionMatrixParam,
                                             Point **&coordinateMeshParam,
                                             TriangleRightPart *&rightPartParam,
                                             SystemParameters &systemParameters);
void createGlobalPressureMatrixHLinear(double **&matrixPressure, TriangleContributionMatrix *&contributionMatrix,
                                       double *&rightPartParam, TriangleRightPart *&localRightPartsParam, int n);
void addBorderConditionsHLinear(double **&matrixResult, double *&rightPartParam, int n,
                                int MATRIX_PRESSURE_SIZE, int OTHER_BORDER, int DOWN_BORDER);

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