#include "FEMTriangles.hpp"

void createLocalContributionMatrixForHConst(TriangleContributionMatrix localMatrix,
                                            Point pointI,
                                            Point pointJ,
                                            Point pointK)
{
    double valB, valG;
    double lambdaX = 1.0;
    double lambdaY = 1.0;
    double *b = new double[3];
    double *c = new double[3];

    double area = countArea(pointI, pointJ, pointK);

    b[0] = pointJ.getY() - pointK.getY();
    b[1] = pointK.getY() - pointI.getY();
    b[2] = pointI.getY() - pointJ.getY();

    c[0] = pointK.getX() - pointJ.getX();
    c[1] = pointI.getX() - pointK.getX();
    c[2] = pointJ.getX() - pointI.getX();

    for (int i = 0; i < 3; i++)
    {
        b[i] = b[i] / (2.0 * area);
        c[i] = c[i] / (2.0 * area);
    }

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            valB = lambdaX * area * b[i] * b[j];
            valG = lambdaY * area * c[i] * c[j];
            localMatrix.setElement(i, j, valB + valG);
        }
    }
}

void createLocalMatrixForEveryElementHConst(TriangleContributionMatrix *&contributionMatrixParam,
                                            Point **&coordinateMeshParam,
                                            int n)
{
    int finiteElementNumber = -1;
    for (int i = 0; i < n - 1; i++)
    {
        for (int j = 0; j < n - 1; j++)
        {
            finiteElementNumber++;
            createLocalContributionMatrixForHConst(contributionMatrixParam[finiteElementNumber],
                                                   coordinateMeshParam[i][j],
                                                   coordinateMeshParam[i + 1][j],
                                                   coordinateMeshParam[i][j + 1]);
            finiteElementNumber++;
            createLocalContributionMatrixForHConst(contributionMatrixParam[finiteElementNumber],
                                                   coordinateMeshParam[i + 1][j + 1],
                                                   coordinateMeshParam[i][j + 1],
                                                   coordinateMeshParam[i + 1][j]);
        }
    }
}

void createGlobalPressureMatrixHConst(double **&matrixPressure,
                                      TriangleContributionMatrix *&contributionMatrix,
                                      int n)
{
    int *globalNodeNumbersIJK = new int[3];
    int finiteElementNumber = 0;

    for (int i = 0; i < n - 1; i++)
    {
        for (int j = 0; j < n - 1; j++)
        {
            //for top triangle
            /*i*/ globalNodeNumbersIJK[0] = i * n + j;
            /*j*/ globalNodeNumbersIJK[1] = globalNodeNumbersIJK[0] + n;
            /*k*/ globalNodeNumbersIJK[2] = globalNodeNumbersIJK[0] + 1;

            for (int iterator1 = 0; iterator1 < 3; iterator1++)
                for (int iterator2 = 0; iterator2 < 3; iterator2++)
                {
                    matrixPressure[globalNodeNumbersIJK[iterator1]][globalNodeNumbersIJK[iterator2]] +=
                        contributionMatrix[finiteElementNumber].matrix[iterator1][iterator2];
                }

            finiteElementNumber++;

            //for bottom triangle
            globalNodeNumbersIJK[0] = globalNodeNumbersIJK[1] + 1;
            std::swap(globalNodeNumbersIJK[1], globalNodeNumbersIJK[2]);

            for (int iterator1 = 0; iterator1 < 3; iterator1++)
                for (int iterator2 = 0; iterator2 < 3; iterator2++)
                {
                    matrixPressure[globalNodeNumbersIJK[iterator1]][globalNodeNumbersIJK[iterator2]] +=
                        contributionMatrix[finiteElementNumber].matrix[iterator1][iterator2];
                }
            finiteElementNumber++;
        }
    }
}

void solveWithHConst(TriangleContributionMatrix *&contributionMatrix,
                     Point **&coordinateMesh,
                     double **&matrixPressure,
                     double *&rightPart,
                     SystemParameters &systemParameters)
{
    const int MATRIX_PRESSURE_SIZE = systemParameters.n * systemParameters.n;
    const int MATRIX_CONTRIBUTION_SIZE = (systemParameters.n - 1) * (systemParameters.n - 1) * 2;

    initMatrix(matrixPressure, MATRIX_PRESSURE_SIZE, MATRIX_PRESSURE_SIZE);
    for (int i = 0; i < MATRIX_PRESSURE_SIZE; i++)
        for (int j = 0; j < MATRIX_PRESSURE_SIZE; j++)
            matrixPressure[i][j] = 0.0;
    initMesh(coordinateMesh, systemParameters);

    initVector(rightPart, MATRIX_PRESSURE_SIZE);
    for (int i = 0; i < MATRIX_PRESSURE_SIZE; i++)
        rightPart[i] = 0.0;

    initContributionMatrix(contributionMatrix, MATRIX_CONTRIBUTION_SIZE);

    createLocalMatrixForEveryElementHConst(contributionMatrix, coordinateMesh, systemParameters.n);
    createGlobalPressureMatrixHConst(matrixPressure, contributionMatrix, systemParameters.n);

    addBorderConditionsOnPressureValues(matrixPressure, rightPart, systemParameters.n, MATRIX_PRESSURE_SIZE,
                                        systemParameters.LOW_BORDER, systemParameters.HIGH_BORDER);

    outputPressureMatrix(matrixPressure, MATRIX_PRESSURE_SIZE);
}

int solveWithHConstBCLR(TriangleContributionMatrix *&contributionMatrix,
                        Point **&coordinateMesh,
                        double **&matrixPressure,
                        SystemParameters &systemParameters)
{
    const int MATRIX_PRESSURE_SIZE = systemParameters.n * systemParameters.n;
    const int MATRIX_CONTRIBUTION_SIZE = (systemParameters.n - 1) * (systemParameters.n - 1) * 2;

    initMatrix(matrixPressure, MATRIX_PRESSURE_SIZE, MATRIX_PRESSURE_SIZE);
    for (int i = 0; i < MATRIX_PRESSURE_SIZE; i++)
        for (int j = 0; j < MATRIX_PRESSURE_SIZE; j++)
            matrixPressure[i][j] = 0.0;

    initMesh(coordinateMesh, systemParameters);

    initContributionMatrix(contributionMatrix, MATRIX_CONTRIBUTION_SIZE);

    createLocalMatrixForEveryElementHConst(contributionMatrix, coordinateMesh, systemParameters.n);
    createGlobalPressureMatrixHConst(matrixPressure, contributionMatrix, systemParameters.n);

    addBorderConditionsToLeftAndRight(matrixPressure, systemParameters.n, 0.0, MATRIX_PRESSURE_SIZE,
                                      systemParameters.LOW_BORDER, systemParameters.HIGH_BORDER);

    outputPressureMatrix(matrixPressure, MATRIX_PRESSURE_SIZE);
    return 0;
}
