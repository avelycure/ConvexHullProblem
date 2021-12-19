#include "FEMTriangles.hpp"

void solveWithTrianglesSecondOrder(TriangleContributionMatrixSecondOrder *&contributionMatrix,
                                   TriangleRightPartSecondOrder *&localRigthParts,
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

    contributionMatrix = new TriangleContributionMatrixSecondOrder[MATRIX_CONTRIBUTION_SIZE];

    localRigthParts = new TriangleRightPartSecondOrder[MATRIX_CONTRIBUTION_SIZE];

    initVector(rightPart, MATRIX_PRESSURE_SIZE);

    createLocalMatrixForEveryElementQuadraticTriangles(contributionMatrix,
                                                       coordinateMesh,
                                                       localRigthParts,
                                                       systemParameters);

    createGlobalPressureMatrixQuadraticTriangles(matrixPressure,
                                                 contributionMatrix,
                                                 rightPart,
                                                 localRigthParts,
                                                 systemParameters.n);
                                                 
    addBorderConditionsQuadraticTriangles(matrixPressure,
                                          rightPart,
                                          systemParameters.n,
                                          MATRIX_PRESSURE_SIZE,
                                          systemParameters.LOW_BORDER,
                                          systemParameters.HIGH_BORDER);

    outputPressureMatrix(matrixPressure, MATRIX_PRESSURE_SIZE);
}

void createLocalMatrixForEveryElementQuadraticTriangles(TriangleContributionMatrixSecondOrder *&contributionMatrixParam,
                                                        Point **&coordinateMeshParam,
                                                        TriangleRightPartSecondOrder *&rightPartParam,
                                                        SystemParameters &systemParameters)
{
    double k = systemParameters.k;
    double hMin = systemParameters.hMin;
    int n = systemParameters.n;

    int finiteElementNumber = -1;
    for (int i = 0; i < n - 2; i += 2)
    {
        for (int j = 0; j < n - 2; j += 2)
        {
            finiteElementNumber++;
            createLocalContributionMatrixForQuardaticTriangleTop(contributionMatrixParam[finiteElementNumber],
                                                                 coordinateMeshParam[i][j],
                                                                 coordinateMeshParam[i + 1][j],
                                                                 coordinateMeshParam[i + 2][j],
                                                                 coordinateMeshParam[i + 1][j + 1],
                                                                 coordinateMeshParam[i][j + 2],
                                                                 coordinateMeshParam[i][j + 1],
                                                                 rightPartParam[finiteElementNumber],
                                                                 systemParameters);
            finiteElementNumber++;
            createLocalContributionMatrixForQuardaticTriangleBottom(contributionMatrixParam[finiteElementNumber],
                                                                    coordinateMeshParam[i + 2][j + 2],
                                                                    coordinateMeshParam[i + 1][j + 2],
                                                                    coordinateMeshParam[i][j + 2],
                                                                    coordinateMeshParam[i + 1][j + 1],
                                                                    coordinateMeshParam[i + 2][j],
                                                                    coordinateMeshParam[i + 2][j + 1],
                                                                    rightPartParam[finiteElementNumber],
                                                                    systemParameters);
        }
    }
}

void createLocalContributionMatrixForQuardaticTriangleTop(TriangleContributionMatrixSecondOrder localMatrix,
                                                          Point pointI,
                                                          Point pointJ,
                                                          Point pointK,
                                                          Point pointL,
                                                          Point pointM,
                                                          Point pointN,
                                                          TriangleRightPartSecondOrder localRightPart,
                                                          SystemParameters &systemParameters)
{
}

void createLocalContributionMatrixForQuardaticTriangleBottom(TriangleContributionMatrixSecondOrder localMatrix,
                                                             Point pointI,
                                                             Point pointJ,
                                                             Point pointK,
                                                             Point pointL,
                                                             Point pointM,
                                                             Point pointN,
                                                             TriangleRightPartSecondOrder localRightPart,
                                                             SystemParameters &systemParameters)
{
}

void createGlobalPressureMatrixQuadraticTriangles(double **&matrixPressure,
                                                  TriangleContributionMatrixSecondOrder *&contributionMatrix,
                                                  double *&rightPartParam,
                                                  TriangleRightPartSecondOrder *&localRightPartsParam,
                                                  int n)
{
    int *globalNodeNumbersIJK = new int[6];
    int finiteElementNumber = 0;

    for (int i = 0; i < n - 2; i += 2)
        for (int j = 0; j < n - 2; j += 2)
        {
            //for top triangle
            /*i*/ globalNodeNumbersIJK[0] = i * n + j;
            /*j*/ globalNodeNumbersIJK[1] = globalNodeNumbersIJK[0] + n;
            /*k*/ globalNodeNumbersIJK[2] = globalNodeNumbersIJK[1] + n;
            /*l*/ globalNodeNumbersIJK[3] = globalNodeNumbersIJK[1] + 1;
            /*m*/ globalNodeNumbersIJK[4] = globalNodeNumbersIJK[0] + 2;
            /*n*/ globalNodeNumbersIJK[5] = globalNodeNumbersIJK[0] + 1;

            for (int i1 = 0; i1 < 6; i1++)
                for (int i2 = 0; i2 < 6; i2++)
                    matrixPressure[globalNodeNumbersIJK[i1]][globalNodeNumbersIJK[i2]] +=
                        contributionMatrix[finiteElementNumber].matrix[i1][i2];

            for (int i1 = 0; i1 < 6; i1++)
                rightPartParam[globalNodeNumbersIJK[i1]] += localRightPartsParam[finiteElementNumber].getElement(i1);

            finiteElementNumber++;

            //for bottom triangle
            /*i*/ globalNodeNumbersIJK[0] = globalNodeNumbersIJK[2] + 2;
            /*j*/ globalNodeNumbersIJK[1] = globalNodeNumbersIJK[0] - n;
            /*k*/ globalNodeNumbersIJK[2] = globalNodeNumbersIJK[1] - n;
            /*l*/ globalNodeNumbersIJK[3] = globalNodeNumbersIJK[1] - 1;
            /*m*/ globalNodeNumbersIJK[4] = globalNodeNumbersIJK[0] - 2;
            /*n*/ globalNodeNumbersIJK[5] = globalNodeNumbersIJK[0] - 1;

            for (int i1 = 0; i1 < 6; i1++)
                for (int i2 = 0; i2 < 6; i2++)
                    matrixPressure[globalNodeNumbersIJK[i1]][globalNodeNumbersIJK[i2]] +=
                        contributionMatrix[finiteElementNumber].matrix[i1][i2];

            for (int i1 = 0; i1 < 6; i1++)
                rightPartParam[globalNodeNumbersIJK[i1]] += localRightPartsParam[finiteElementNumber].getElement(i1);

            finiteElementNumber++;
        }
}

void addBorderConditionsQuadraticTriangles(double **&matrixResult,
                                           double *&rightPartParam,
                                           int n,
                                           int MATRIX_PRESSURE_SIZE,
                                           int OTHER_BORDER,
                                           int DOWN_BORDER)
{
    //0 row
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < MATRIX_PRESSURE_SIZE; j++)
        {
            matrixResult[i][j] = 0.0;
        }
        matrixResult[i][i] = 1.0;
        rightPartParam[i] = DOWN_BORDER;
    }

    //left
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < MATRIX_PRESSURE_SIZE; j++)
        {
            matrixResult[i * n][j] = 0.0;
        }
        matrixResult[i * n][i * n] = 1.0;
        rightPartParam[i * n] = OTHER_BORDER;
    }

    //right
    for (int i = 1; i <= n; i++)
    {
        for (int j = 0; j < MATRIX_PRESSURE_SIZE; j++)
        {
            matrixResult[i * n - 1][j] = 0.0;
        }
        matrixResult[i * n - 1][i * n - 1] = 1.0;
        rightPartParam[i * n - 1] = OTHER_BORDER;
    }

    //n row
    for (int i = MATRIX_PRESSURE_SIZE - n; i < MATRIX_PRESSURE_SIZE; i++)
    {
        for (int j = 0; j < MATRIX_PRESSURE_SIZE; j++)
        {
            matrixResult[i][j] = 0.0;
        }
        matrixResult[i][i] = 1.0;
        rightPartParam[i] = OTHER_BORDER;
    }

    //displayMatrix(matrixResult, MATRIX_PRESSURE_SIZE, MATRIX_PRESSURE_SIZE);

    std::fstream myFile;
    myFile.open("data/fem_output/rightPart.txt", std::fstream::out);
    for (int i = 0; i < MATRIX_PRESSURE_SIZE; i++)
        myFile << rightPartParam[i] << std::endl;
}