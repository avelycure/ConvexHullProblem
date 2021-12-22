#include "InitFuncs.hpp"

void initMatrix(double **&matrix, int row, int column)
{
    matrix = new double *[row];
    for (int i = 0; i < row; i++)
        matrix[i] = new double[column];
}

void displayMatrix(double **matrix, int row, int column)
{
    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < column; j++)
            std::cout << matrix[i][j] << " ";
        std::cout << std::endl;
    }
}

/**
 * Show X and Y coordinates of the nodes of the mesh
 * */
void displayMesh(Point **coordinateMesh, int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
            std::cout << "(" << coordinateMesh[i][j].getX() << "," << coordinateMesh[i][j].getY() << ") ";
        std::cout << std::endl;
    }
}

void displayVector(double *mVector, int n)
{
    for (int i = 0; i < n; i++)
        std::cout << mVector[i] << " ";
    std::cout << std::endl;
}

void displayAllLocalMatrixes(TriangleContributionMatrix *&contributionMatrixParam, int n)
{
    for (int i = 0; i < n; i++)
    {
        std::cout << "matrix[" << i << "]" << std::endl;
        displayMatrix(contributionMatrixParam[i].matrix, 3, 3);
        std::cout << std::endl;
    }
}

void outputPressureMatrix(double **matrixPressure, int MATRIX_PRESSURE_SIZE)
{
    std::fstream myFile;

    myFile.open("data/fem_output/matrixPressure.txt", std::fstream::out);
    for (int i = 0; i < MATRIX_PRESSURE_SIZE; i++)
    {
        for (int j = 0; j < MATRIX_PRESSURE_SIZE; j++)
            myFile << matrixPressure[i][j] << " ";
        myFile << std::endl;
    }
}

void initVector(double *&p, int n)
{
    p = new double[n];
}

void initContributionMatrix(TriangleContributionMatrix *&contributionMatrix, int MATRIX_CONTRIBUTION_SIZE)
{
    contributionMatrix = new TriangleContributionMatrix[MATRIX_CONTRIBUTION_SIZE];
}

void initRightPart(TriangleRightPart *&localRigthParts, int MATRIX_CONTRIBUTION_SIZE)
{
    localRigthParts = new TriangleRightPart[MATRIX_CONTRIBUTION_SIZE];
}

void initMesh(Point **&coordinateMesh, SystemParameters &systemParameters)
{
    int n = systemParameters.n;
    double xOrigin = systemParameters.xOrigin;
    double yOrigin = systemParameters.yOrigin;
    double h = systemParameters.borderLength / (n - 1);
    coordinateMesh = new Point *[n];
    for (int i = 0; i < n; i++)
        coordinateMesh[i] = new Point[n];

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
        {
            coordinateMesh[i][j].setX(xOrigin + i * h);
            coordinateMesh[i][j].setY(yOrigin + j * h);
        }
}

void readSystemParameters(SystemParameters &systemParameters, std::string &method)
{
    nlohmann::json j;
    std::fstream fileInputSystem;
    fileInputSystem.open(FILE_SYSTEM_NAME);
    fileInputSystem >> j;
    int k = j["system"];
    fileInputSystem.close();

    nlohmann::json jsonObject;
    std::fstream fileInput;
    fileInput.open(FILE_PARAMETERS_NAME);
    fileInput >> jsonObject;

    method = jsonObject["configs"][k]["method"];
    systemParameters.n = jsonObject["configs"][k]["n"];
    systemParameters.L = jsonObject["configs"][k]["L"];
    systemParameters.U = jsonObject["configs"][k]["U"];
    systemParameters.k = jsonObject["configs"][k]["k"];
    systemParameters.mu = jsonObject["configs"][k]["mu"];
    systemParameters.Hn = jsonObject["configs"][k]["Hn"]; // для обезразмеривания
    systemParameters.pMin = jsonObject["configs"][k]["pMin"];
    systemParameters.hMin = jsonObject["configs"][k]["hMin"]; //высота впускного зазора
    systemParameters.xOrigin = jsonObject["configs"][k]["xOrigin"];
    systemParameters.yOrigin = jsonObject["configs"][k]["yOrigin"];
    systemParameters.LOW_BORDER = jsonObject["configs"][k]["lowBorder"];
    systemParameters.HIGH_BORDER = jsonObject["configs"][k]["highBorder"];
    systemParameters.borderLength = jsonObject["configs"][k]["borderLength"];

    fileInput.close();
}
