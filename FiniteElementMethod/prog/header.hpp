#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <cmath>
#include <sstream>
#include <algorithm>
#include "Point.h"
#include "ContributionMatrix.h"
#include "SystemParameters.hpp"
#include "single_include/nlohmann/json.hpp"
using namespace std;

const string FILE_PARAMETERS_NAME = "systemParameters.json";
const string FILE_SYSTEM_NAME = "systemNum.json";

/**
 * Name of file for matrix of pressures in node
 * */
const string FILE_OUTPUT_MATRIX = "data/pressureMatrix.txt";

/**
 * Name of file for right part
 * */
const string FILE_VECTOR_RIGHT_PART_OUTPUT = "data/rightPart.txt";

/**
 * Names of methods
 * */
const string H_CONST = "hConst";
const string H_LINEAR = "hLinear";

/**
 * 
 * 
 * 
 * Common functions
 * 
 * 
 * 
 * */

void initMatrix(double **&matrix, int row, int column);
void displayMatrix(double **matrix, int row, int column);
void displayMesh(Point **coordinateMesh, int n);
void displayAllLocalMatrixes(ContributionMatrix *&ContributionMatrixParam, int n);
void outputPressureMatrix(double **matrixPressure, int MATRIX_PRESSURE_SIZE);
void displayVector(double *mVector, int n);
void initVector(double *&p, int n);
void initMesh(Point **&coordinateMesh, SystemPatemeters &systemParameters);
void initRightPart(RightPart *&localRigthParts, int MATRIX_CONTRIBUTION_SIZE);
void initContributionMatrix(ContributionMatrix *&contributionMatrix, int MATRIX_CONTRIBUTION_SIZE);
void readSystemParameters(SystemPatemeters &systemParameters, string &method);
void dimensionlessSystemParameters(SystemPatemeters &systemParameters, string &method);

/**
 * 
 * 
 * 
 * Finite elements method functions
 * 
 * 
 * 
 * */

/**
 * Constant height
 * */
void solveWithHConst(ContributionMatrix *&contributionMatrix,
                    Point **&coordinateMesh,
                    double **&matrixPressure,
                    SystemPatemeters &systemParameters);
double countArea(Point pointI, Point pointJ, Point pointK);
void createLocalContributionMatrixForHConst(ContributionMatrix localMatrix, Point pointI, Point pointJ, Point pointK);
void createLocalMatrixForEveryElementHConst(ContributionMatrix *&contributionMatrixParam, Point **&coordinateMeshParam, int n);
void addBorderConditionsHConst(double **&matrixResult, int n, int MATRIX_PRESSURE_SIZE, int OTHER_BORDER, int DOWN_BORDER);
void createGlobalPressureMatrixHConst(double **&matrixPressure, ContributionMatrix *&contributionMatrix, int n);

/**
 * Linear changing height
 * */
void solveWithHLinear(ContributionMatrix *&contributionMatrix,
                     RightPart *&localRigthParts,
                     Point **&coordinateMesh,
                     double **&matrixPressure,
                     double *&rightPart,
                     SystemPatemeters &systemParameters);
void createLocalContributionMatrixForHLinearTop(ContributionMatrix localMatrix,
                                               Point pointI, Point pointJ, Point pointK,
                                               RightPart localRightPart,SystemPatemeters &systemParameters);
void createLocalContributionMatrixForHLinearBottom(ContributionMatrix localMatrix,
                                                  Point pointI, Point pointJ, Point pointK,
                                                  RightPart localRightPart,SystemPatemeters &systemParameters);
void createLocalMatrixForEveryElementHLinear(ContributionMatrix *&contributionMatrixParam,
                                             Point **&coordinateMeshParam,
                                             RightPart *&rightPartParam,
                                             SystemPatemeters &systemParameters);
void createGlobalPressureMatrixHLinear(double **&matrixPressure, ContributionMatrix *&contributionMatrix,
                                       double *&rightPartParam, RightPart *&localRightPartsParam, int n);
void addBorderConditionsHLinear(double **&matrixResult, double *&rightPartParam, int n,
                                int MATRIX_PRESSURE_SIZE, int OTHER_BORDER, int DOWN_BORDER);

/**
 * 
 * 
 * 
 * Gauss functions
 * 
 * 
 * 
 * */
int solveEquation(const int size);
int allocateMemory(double **&A, double *&B, double *&X1, const int &n);
int readData(const string fileNameMatrix, const string fileNameVector, double **&matrixA, double *&vectorB, const int &n);
void writeVector(string fileNameOutput, double *&vector, const int &n);
bool gaussMethod(double **&A, double *&B, double *&X, const int &size);
bool matrixIsPrepared(double **&A, double *&B, const int &i, vector<tuple<int, int>> &permutations, const int &size);
void diagonalizeEquation(double **&A, double *&B, double *&X, const int &size, vector<tuple<int, int>> &permutations);
tuple<int, int> searchMax(double **&A, const int &currentMinor, const int &size);
void swapRows(double **&A, double *&B, const int &i1, const int &i2);
void swapColomns(double **&A, const int &j1, const int &j2, const int &size);
bool isDegenerate(double **&A, const int &i, const int &size);
void freeMemory(double **&A, double *&B, double *&X1, const int &n);