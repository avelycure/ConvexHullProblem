#include "SecondOrderRectangleFE.hpp"

void solveWithSecondOrderRectangleFE(SecondOrderRectangleContributionMatrix *&contributionMatrix,
                                     SecondOrderRectangleRightPart *&localRigthParts,
                                     Point **&coordinateMesh,
                                     double **&matrixPressure,
                                     double *&rightPart,
                                     SystemParameters &systemParameters)
{
    int numberOfFE = ((systemParameters.n - 1) / 2) * ((systemParameters.n - 1) / 2);

    //substract number of points in center
    const int MATRIX_PRESSURE_SIZE = systemParameters.n * systemParameters.n - numberOfFE;
    const int MATRIX_CONTRIBUTION_SIZE = ((systemParameters.n - 1) / 2) * ((systemParameters.n - 1) / 2);

    initMatrix(matrixPressure, MATRIX_PRESSURE_SIZE, MATRIX_PRESSURE_SIZE);
    for (int i = 0; i < MATRIX_PRESSURE_SIZE; i++)
        for (int j = 0; j < MATRIX_PRESSURE_SIZE; j++)
            matrixPressure[i][j] = 0.0;
    initMesh(coordinateMesh, systemParameters);

    contributionMatrix = new SecondOrderRectangleContributionMatrix[MATRIX_CONTRIBUTION_SIZE];

    localRigthParts = new SecondOrderRectangleRightPart[MATRIX_CONTRIBUTION_SIZE];

    rightPart = new double[MATRIX_PRESSURE_SIZE];
    for (int i = 0; i < MATRIX_PRESSURE_SIZE; i++)
        rightPart[i] = 0.0;
}

void createLocalMatrixForEveryRectangleElementSecondOrder(SecondOrderRectangleContributionMatrix *&contributionMatrixParam,
                                                          Point **&coordinateMeshParam,
                                                          SecondOrderRectangleRightPart *&rightPartParam,
                                                          SystemParameters &systemParameters)
{
    int n = systemParameters.n;
    int element = 0;

    for (int i = 0; i < n - 2; i += 2)
        for (int j = 0; j < n - 2; j += 2)
        {
            createLocalContributionMatrixForRectangleElementSecondOrder(contributionMatrixParam[element],
                                                                        coordinateMeshParam[i][j],
                                                                        coordinateMeshParam[i + 1][j],
                                                                        coordinateMeshParam[i + 2][j],
                                                                        coordinateMeshParam[i + 2][j + 1],
                                                                        coordinateMeshParam[i + 2][j + 2],
                                                                        coordinateMeshParam[i + 1][j + 2],
                                                                        coordinateMeshParam[i][j + 2],
                                                                        coordinateMeshParam[i][j + 1],
                                                                        rightPartParam[element],
                                                                        systemParameters);
            element++;
        }
}

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
                                                                 SystemParameters &systemParameters)
{
    double hMin = systemParameters.hMin;
    double k = systemParameters.k;

    //coefficients of polynom
    double A1 = pow(hMin, 3.0);
    double A2 = 3.0 * pow(hMin, 2.0) * k;
    double A3 = 3.0 * hMin * pow(k, 2.0);
    double A4 = pow(k, 3.0);

    double s1 = pointR.getY() - pointI.getY();
    double s2 = (1.0 / 2.0) * (pow(pointR.getY(), 2.0) - pow(pointI.getY(), 2.0));
    double s3 = (1.0 / 3.0) * (pow(pointR.getY(), 3.0) - pow(pointI.getY(), 3.0));
    double s4 = (1.0 / 4.0) * (pow(pointR.getY(), 4.0) - pow(pointI.getY(), 4.0));
    double s5 = (1.0 / 5.0) * (pow(pointR.getY(), 5.0) - pow(pointI.getY(), 5.0));
    double s6 = (1.0 / 6.0) * (pow(pointR.getY(), 6.0) - pow(pointI.getY(), 6.0));
    double s7 = (1.0 / 7.0) * (pow(pointR.getY(), 7.0) - pow(pointI.getY(), 7.0));
    double s8 = (1.0 / 8.0) * (pow(pointR.getY(), 8.0) - pow(pointI.getY(), 8.0));

    double z1 = pointK.getX() - pointI.getX();
    double z2 = (1.0 / 2.0) * (pow(pointK.getX(), 2.0) - pow(pointI.getX(), 2.0));
    double z3 = (1.0 / 3.0) * (pow(pointK.getX(), 3.0) - pow(pointI.getX(), 3.0));
    double z4 = (1.0 / 4.0) * (pow(pointK.getX(), 4.0) - pow(pointI.getX(), 4.0));
    double z5 = (1.0 / 5.0) * (pow(pointK.getX(), 5.0) - pow(pointI.getX(), 5.0));

    //clR1
    double c1 = A1 * s1 * z1 + A2 * s2 * z1 + A3 * s3 * z1 + A4 * s4 * z1;

    //blR6
    double c2 = A1 * s1 * z1 + A2 * s2 * z1 + A3 * s3 * z1 + A4 * s4 * z1;

    //blR2
    double c3 = A1 * s2 * z1 + A2 * s3 * z1 + A3 * s4 * z1 + A4 * s5 * z1;

    //dlR6
    double c4 = A1 * s2 * z1 + A2 * s3 * z1 + A3 * s4 * z1 + A4 * s5 * z1;

    //flR1
    double c5 = 2 * A1 * s2 * z1 + 2 * A2 * s3 * z1 + 2 * A3 * s4 * z1 + 2 * A4 * s5 * z1;

    //clR1
    double c6 = 2 * A1 * s2 * z1 + 2 * A2 * s3 * z1 + 2 * A3 * s4 * z1 + 2 * A4 * s5 * z1;

    //blR5
    double c7 = A1 * s3 * z1 + A2 * s4 * z1 + A3 * s5 * z1 + A4 * s6 * z1;

    //tlR6
    double c8 = A1 * s3 * z1 + A2 * s4 * z1 + A3 * s5 * z1 + A4 * s6 * z1;

    //flR3
    double c9 = 4 * A1 * s3 * z1 + 4 * A2 * s4 * z1 + 4 * A3 * s5 * z1 + 4 * A4 * s6 * z1;

    //dlR1
    double c10 = A1 * s1 * z2 + A2 * s2 * z2 + A3 * s3 * z2 + A4 * s4 * z2;

    //clR2
    double c11 = A1 * s1 * z2 + A2 * s2 * z2 + A3 * s3 * z2 + A4 * s4 * z2;

    //elR6
    double c12 = 2 * A1 * s2 * z2 + 2 * A2 * s3 * z2 + 2 * A3 * s4 * z2 + 2 * A4 * s5 * z2;

    //flR2
    double c13 = 2 * A1 * s2 * z2 + 2 * A2 * s3 * z2 + 2 * A3 * s4 * z2 + 2 * A4 * s5 * z2;

    //dlR3
    double c14 = 2 * A1 * s2 * z2 + 2 * A2 * s3 * z2 + 2 * A3 * s4 * z2 + 2 * A4 * s5 * z2;

    //blR4
    double c15 = 2 * A1 * s2 * z2 + 2 * A2 * s3 * z2 + 2 * A3 * s4 * z2 + 2 * A4 * s5 * z2;

    //clR5
    double c16 = 2 * A1 * s2 * z2 + 2 * A2 * s3 * z2 + 2 * A3 * s4 * z2 + 2 * A4 * s5 * z2;

    //glR6
    double c17 = 2 * A1 * s2 * z2 + 2 * A2 * s3 * z2 + 2 * A3 * s4 * z2 + 2 * A4 * s5 * z2;

    //dlR7
    double c18 = 2 * A1 * s2 * z2 + 2 * A2 * s3 * z2 + 2 * A3 * s4 * z2 + 2 * A4 * s5 * z2;

    //tlR1
    double c19 = 2 * A1 * s2 * z2 + 2 * A2 * s3 * z2 + 2 * A3 * s4 * z2 + 2 * A4 * s5 * z2;

    //elR5
    double c20 = 2 * A1 * s3 * z2 + 2 * A2 * s4 * z2 + 2 * A3 * s5 * z2 + 2 * A4 * s6 * z2;

    //tlR7
    double c21 = 2 * A1 * s3 * z2 + 2 * A2 * s4 * z2 + 2 * A3 * s5 * z2 + 2 * A4 * s6 * z2;

    //flR5
    double c22 = 4 * A1 * s3 * z2 + 4 * A2 * s4 * z2 + 4 * A3 * s5 * z2 + 4 * A4 * s6 * z2;

    //tlR3
    double c23 = 4 * A1 * s3 * z2 + 4 * A2 * s4 * z2 + 4 * A3 * s5 * z2 + 4 * A4 * s6 * z2;

    //glR1
    double c24 = A1 * s1 * z3 + A2 * s2 * z3 + A3 * s3 * z3 + A4 * s4 * z3;

    //clR4
    double c25 = A1 * s1 * z3 + A2 * s2 * z3 + A3 * s3 * z3 + A4 * s4 * z3;

    //dlR2
    double c26 = A1 * s3 * z1 + A2 * s4 * z1 + A3 * s5 * z1 + A4 * s6 * z1 +
                 A1 * s1 * z3 + A2 * s2 * z3 + A3 * s3 * z3 + A4 * s4 * z3;

    //elR7
    double c27 = 4 * A1 * s1 * z3 + 4 * A2 * s2 * z3 + 4 * A3 * s3 * z3 + 4 * A4 * s4 * z3;

    //flR4
    double c28 = 2 * A1 * s2 * z3 + 2 * A2 * s3 * z3 + 2 * A3 * s4 * z3 + 2 * A4 * s5 * z3;

    //dlR5
    double c29 = A1 * s4 * z1 + A2 * s5 * z1 + A3 * s6 * z1 + A4 * s7 * z1 +
                 2 * A1 * s2 * z3 + 2 * A2 * s3 * z3 + 2 * A3 * s4 * z3 + 2 * A4 * s5 * z3;

    //tlR2
    double c30 = A1 * s4 * z1 + A2 * s5 * z1 + A3 * s6 * z1 + A4 * s7 * z1 +
                 2 * A1 * s2 * z3 + 2 * A2 * s3 * z3 + 2 * A3 * s4 * z3 + 2 * A4 * s5 * z3;

    //elR4
    double c31 = 4 * A1 * s2 * z3 + 4 * A2 * s3 * z3 + 4 * A3 * s4 * z3 + 4 * A4 * s5 * z3;

    //glR7
    double c32 = 4 * A1 * s2 * z3 + 4 * A2 * s3 * z3 + 4 * A3 * s4 * z3 + 4 * A4 * s5 * z3;

    //glR3
    double c33 = 2 * A1 * s2 * z3 + 2 * A2 * s3 * z3 + 2 * A3 * s4 * z3 +
                 2 * A4 * s5 * z3;

    //tlR5
    double c34 = A1 * s5 * z1 + A2 * s6 * z1 + A3 * s7 * z1 + A4 * s8 * z1 +
                 4 * A1 * s3 * z3 + 4 * A2 * s4 * z3 + 4 * A3 * s5 * z3 + 4 * A4 * s6 * z3;

    //glR2
    double c35 = 2 * A1 * s3 * z2 + 2 * A2 * s4 * z2 + 2 * A3 * s5 * z2 + 2 * A4 * s6 * z2 +
                 A1 * s1 * z4 + A2 * s2 * z4 + A3 * s3 * z4 + A4 * s4 * z4;

    //dlR4
    double c36 = 2 * A1 * s3 * z2 + 2 * A2 * s4 * z2 + 2 * A3 * s5 * z2 + 2 * A4 * s6 * z2 +
                 A1 * s1 * z4 + A2 * s2 * z4 + A3 * s3 * z4 + A4 * s4 * z4;

    //glR5
    double c37 = 2 * A1 * s4 * z2 + 2 * A2 * s5 * z2 + 2 * A3 * s6 * z2 + 2 * A4 * s7 * z2 +
                 2 * A1 * s2 * z4 + 2 * A2 * s3 * z4 + 2 * A3 * s4 * z4 + 2 * A4 * s5 * z4;

    //tlR4
    double c38 = 2 * A1 * s4 * z2 + 2 * A2 * s5 * z2 + 2 * A3 * s6 * z2 + 2 * A4 * s7 * z2 +
                 2 * A1 * s2 * z4 + 2 * A2 * s3 * z4 + 2 * A3 * s4 * z4 + 2 * A4 * s5 * z4;

    //glR4
    double c39 = 4 * A1 * s3 * z3 + 4 * A2 * s4 * z3 + 4 * A3 * s5 * z3 + 4 * A4 * s6 * z3 +
                 A1 * s1 * z5 + A2 * s2 * z5 + A3 * s3 * z5 + A4 * s4 * z5;
}