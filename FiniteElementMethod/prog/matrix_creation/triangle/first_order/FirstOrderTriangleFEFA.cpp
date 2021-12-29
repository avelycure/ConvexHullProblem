#include "FirstOrderTriangleFEFA.hpp"

/**
 * Functional approach. These mehtod gives solution only for constant height
 * */

void solveWithFirstOrderTriangleFEConstantHeight(TriangleContributionMatrix *&contributionMatrix,
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

    createLocalMatrixForEveryElementConstantHeight(contributionMatrix, coordinateMesh, systemParameters.n);
    createGlobalPressureMatrixConstantHeight(matrixPressure, contributionMatrix, systemParameters.n);

    addBorderConditionsOnPressureValues(matrixPressure, rightPart, systemParameters.n, MATRIX_PRESSURE_SIZE,
                                        systemParameters.LOW_BORDER, systemParameters.HIGH_BORDER);

    outputPressureMatrix(matrixPressure, MATRIX_PRESSURE_SIZE);
}

void createLocalContributionMatrixConstantHeight(TriangleContributionMatrix localMatrix,
                                            Point pointI,
                                            Point pointJ,
                                            Point pointK)
{
    double valB, valG;
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
        for (int j = 0; j < 3; j++)
        {
            valB = area * b[i] * b[j];
            valG = area * c[i] * c[j];
            localMatrix.setElement(i, j, valB + valG);
        }

    delete[] b;
    delete[] c;
}

void createLocalMatrixForEveryElementConstantHeight(TriangleContributionMatrix *&contributionMatrixParam,
                                            Point **&coordinateMeshParam,
                                            int n)
{
    int finiteElementNumber = -1;
    for (int i = 0; i < n - 1; i++)
        for (int j = 0; j < n - 1; j++)
        {
            finiteElementNumber++;
            createLocalContributionMatrixConstantHeight(contributionMatrixParam[finiteElementNumber],
                                                   coordinateMeshParam[i][j],
                                                   coordinateMeshParam[i + 1][j],
                                                   coordinateMeshParam[i][j + 1]);
            finiteElementNumber++;
            createLocalContributionMatrixConstantHeight(contributionMatrixParam[finiteElementNumber],
                                                   coordinateMeshParam[i + 1][j + 1],
                                                   coordinateMeshParam[i][j + 1],
                                                   coordinateMeshParam[i + 1][j]);
        }
}

void createGlobalPressureMatrixConstantHeight(double **&matrixPressure,
                                      TriangleContributionMatrix *&contributionMatrix,
                                      int n)
{
    int *globalNodeNumbersIJK = new int[3];
    int finiteElementNumber = 0;

    for (int i = 0; i < n - 1; i++)
        for (int j = 0; j < n - 1; j++)
        {
            //for top triangle
            /*i*/ globalNodeNumbersIJK[0] = i * n + j;
            /*j*/ globalNodeNumbersIJK[1] = globalNodeNumbersIJK[0] + n;
            /*k*/ globalNodeNumbersIJK[2] = globalNodeNumbersIJK[0] + 1;

            for (int i1 = 0; i1 < 3; i1++)
                for (int i2 = 0; i2 < 3; i2++)
                    matrixPressure[globalNodeNumbersIJK[i1]][globalNodeNumbersIJK[i2]] +=
                        contributionMatrix[finiteElementNumber].matrix[i1][i2];

            finiteElementNumber++;

            //for bottom triangle
            globalNodeNumbersIJK[0] = globalNodeNumbersIJK[1] + 1;
            std::swap(globalNodeNumbersIJK[1], globalNodeNumbersIJK[2]);

            for (int i1 = 0; i1 < 3; i1++)
                for (int i2 = 0; i2 < 3; i2++)
                    matrixPressure[globalNodeNumbersIJK[i1]][globalNodeNumbersIJK[i2]] +=
                        contributionMatrix[finiteElementNumber].matrix[i1][i2];
                
            finiteElementNumber++;
        }

    delete[] globalNodeNumbersIJK;
}