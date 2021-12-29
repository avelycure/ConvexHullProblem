#include "FirstOrderTriangleFEGA.hpp"
/**
  * Galerkin approach
  * */

void solveWithFirstOrderTriangleFE(TriangleContributionMatrix *&contributionMatrix,
                                   TriangleRightPart *&localRigthParts,
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

    initContributionMatrix(contributionMatrix, MATRIX_CONTRIBUTION_SIZE);
    initRightPart(localRigthParts, MATRIX_CONTRIBUTION_SIZE);
    initVector(rightPart, MATRIX_PRESSURE_SIZE);

    createLocalMatrixes(contributionMatrix, coordinateMesh, localRigthParts, systemParameters);

    createGlobalPressureMatrix(matrixPressure, contributionMatrix, rightPart, localRigthParts, systemParameters.n);
    addBorderConditionsOnPressureValues(matrixPressure, rightPart, systemParameters.n, MATRIX_PRESSURE_SIZE,
                                        systemParameters.LOW_BORDER, systemParameters.HIGH_BORDER);

    outputPressureMatrix(matrixPressure, MATRIX_PRESSURE_SIZE);
}

void createLocalContributionMatrix(TriangleContributionMatrix localMatrix,
                                   Point pointI,
                                   Point pointJ,
                                   Point pointK,
                                   TriangleRightPart localRightPart,
                                   SystemParameters &systemParameters)
{
    // coefficients of line z = k1 * x + k2
    double k1 = (pointK.getX() - pointJ.getX()) / (pointK.getY() - pointJ.getY());
    double k2 = pointJ.getX() - pointJ.getY() * (pointK.getX() - pointJ.getX()) / (pointK.getY() - pointJ.getY());

    double hMin = systemParameters.hMin;
    double k = systemParameters.k;

    //coefficients of polynom
    //безразмерное к
    double A1 = pow(hMin, 3.0);
    double A2 = 3.0 * pow(hMin, 2.0) * k;
    double A3 = 3.0 * hMin * pow(k, 2.0);
    double A4 = pow(k, 3.0);

    //
    double s1 = pointK.getY() - pointI.getY();
    double s2 = (1.0 / 2.0) * (pow(pointK.getY(), 2.0) - pow(pointI.getY(), 2.0));
    double s3 = (1.0 / 3.0) * (pow(pointK.getY(), 3.0) - pow(pointI.getY(), 3.0));
    double s4 = (1.0 / 4.0) * (pow(pointK.getY(), 4.0) - pow(pointI.getY(), 4.0));
    double s5 = (1.0 / 5.0) * (pow(pointK.getY(), 5.0) - pow(pointI.getY(), 5.0));

    //blR1
    double c1 = A1 * k1 * s2 + A2 * k2 * s2 + A2 * k1 * s3 + A3 * k2 * s3 +
                A4 * k1 * s5 + A1 * s1 * (k2 - pointI.getX()) - A2 * s2 * pointI.getX() - A3 * s3 * pointI.getX() +
                s4 * (A3 * k1 + A4 * k2 - A4 * pointI.getX());

    //clR2
    double c2 = A1 * k1 * s2 + A2 * k2 * s2 + A2 * k1 * s3 + A3 * k2 * s3 +
                A4 * k1 * s5 + A1 * s1 * (k2 - pointI.getX()) - A2 * s2 * pointI.getX() - A3 * s3 * pointI.getX() +
                s4 * (A3 * k1 + A4 * k2 - A4 * pointI.getX());

    //coefficient R3 6 * mu * k * L * ...
    //тут вроде бы безразмерное к
    double R3 = 6.0 * systemParameters.mu * systemParameters.L * systemParameters.U * k /
                (systemParameters.Hn * systemParameters.Hn * systemParameters.pMin);

    double valR1, valR2;
    double *a = new double[3];
    double *b = new double[3];
    double *c = new double[3];

    a[0] = -(pointJ.getX() / (pointI.getX() - pointJ.getX())) + pointI.getY() / (-pointI.getY() + pointK.getY());
    a[1] = pointI.getX() / (pointI.getX() - pointJ.getX());
    a[2] = pointI.getY() / (pointI.getY() - pointK.getY());

    b[0] = 1.0 / (pointI.getY() - pointK.getY());
    b[1] = 0.0;
    b[2] = 1.0 / (-pointI.getY() + pointK.getY());

    c[0] = 1.0 / (pointI.getX() - pointJ.getX());
    c[1] = 1.0 / (-pointI.getX() + pointJ.getX());
    c[2] = 0.0;

    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
        {
            valR1 = b[i] * b[j] * c1;
            valR2 = c[i] * c[j] * c2;
            localMatrix.setElement(i, j, valR1 + valR2);
        }

    //creating local vector
    for (int i = 0; i < 3; i++)
        localRightPart.setElement(i, -R3 * (a[i] * k1 * s2 + b[i] * k2 * s2 + c[i] * k1 * k2 * s2 + b[i] * k1 * s3 +
                                            0.5 * c[i] * k1 * k1 * s3 + a[i] * s1 * (k2 - pointI.getX()) - b[i] * s2 * pointI.getX() +
                                            0.5 * c[i] * s1 * (k2 * k2 - pointI.getX() * pointI.getX())));

    delete[] a;
    delete[] b;
    delete[] c;
}

void createLocalMatrixes(TriangleContributionMatrix *&contributionMatrixParam,
                         Point **&coordinateMeshParam,
                         TriangleRightPart *&rightPartParam,
                         SystemParameters &systemParameters)
{
    int n = systemParameters.n;
    int finiteElementNumber = -1;
    
    for (int i = 0; i < n - 1; i++)
        for (int j = 0; j < n - 1; j++)
        {
            finiteElementNumber++;
            createLocalContributionMatrix(contributionMatrixParam[finiteElementNumber],
                                          coordinateMeshParam[i][j],
                                          coordinateMeshParam[i + 1][j],
                                          coordinateMeshParam[i][j + 1],
                                          rightPartParam[finiteElementNumber],
                                          systemParameters);
            finiteElementNumber++;
            createLocalContributionMatrix(contributionMatrixParam[finiteElementNumber],
                                          coordinateMeshParam[i + 1][j + 1],
                                          coordinateMeshParam[i][j + 1],
                                          coordinateMeshParam[i + 1][j],
                                          rightPartParam[finiteElementNumber],
                                          systemParameters);
        }
}

void createGlobalPressureMatrix(double **&matrixPressure,
                                TriangleContributionMatrix *&contributionMatrix,
                                double *&rightPartParam,
                                TriangleRightPart *&localRightPartsParam,
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

            for (int i1 = 0; i1 < 3; i1++)
                rightPartParam[globalNodeNumbersIJK[i1]] += localRightPartsParam[finiteElementNumber].getElement(i1);

            finiteElementNumber++;

            //for bottom triangle
            globalNodeNumbersIJK[0] = globalNodeNumbersIJK[1] + 1;
            std::swap(globalNodeNumbersIJK[1], globalNodeNumbersIJK[2]);

            for (int i1 = 0; i1 < 3; i1++)
                for (int i2 = 0; i2 < 3; i2++)
                    matrixPressure[globalNodeNumbersIJK[i1]][globalNodeNumbersIJK[i2]] +=
                        contributionMatrix[finiteElementNumber].matrix[i1][i2];

            for (int i1 = 0; i1 < 3; i1++)
                rightPartParam[globalNodeNumbersIJK[i1]] += localRightPartsParam[finiteElementNumber].getElement(i1);

            finiteElementNumber++;
        }
}