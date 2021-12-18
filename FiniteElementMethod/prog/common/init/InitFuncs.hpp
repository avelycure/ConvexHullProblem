#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "../classes/mesh/Point.hpp"
#include "../classes/contribution_matrix/TriangleContributionMatrix.hpp"
#include "../classes/contribution_matrix/TriangleRightPart.hpp"
#include "../classes/system/SystemParameters.hpp"
#include "../../libs/single_include/nlohmann/json.hpp"

const std::string FILE_PARAMETERS_NAME = "data/fem_input/initial_conditions/systemParameters.json";
const std::string FILE_SYSTEM_NAME = "data/fem_input/initial_conditions/systemNum.json";
const std::string H_CONST = "hConst";
const std::string H_LINEAR = "hLinear";

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
void readSystemParameters(SystemParameters &systemParameters, std::string &method);