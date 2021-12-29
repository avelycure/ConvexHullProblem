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

    createLocalMatrixForEveryRectangleElementSecondOrder(contributionMatrix,
                                                         coordinateMesh,
                                                         localRigthParts,
                                                         systemParameters);

    createGlobalPressureMatrixForRectangleElementSecondOrder(matrixPressure,
                                                             contributionMatrix,
                                                             rightPart,
                                                             localRigthParts,
                                                             systemParameters.n);

    addBorderConditionsForRectangleSecondOrder(matrixPressure,
                                               rightPart,
                                               systemParameters.n,
                                               MATRIX_PRESSURE_SIZE,
                                               systemParameters);

    outputPressureMatrix(matrixPressure, MATRIX_PRESSURE_SIZE);
}

void createLocalMatrixForEveryRectangleElementSecondOrder(SecondOrderRectangleContributionMatrix *&contributionMatrix,
                                                          Point **&coordinateMesh,
                                                          SecondOrderRectangleRightPart *&rightPart,
                                                          SystemParameters &systemParameters)
{
    int n = systemParameters.n;
    int element = 0;

    for (int i = 0; i < n - 2; i += 2)
        for (int j = 0; j < n - 2; j += 2)
        {
            createLocalContributionMatrixForRectangleElementSecondOrder(contributionMatrix[element],
                                                                        coordinateMesh[i][j],
                                                                        coordinateMesh[i + 1][j],
                                                                        coordinateMesh[i + 2][j],
                                                                        coordinateMesh[i + 2][j + 1],
                                                                        coordinateMesh[i + 2][j + 2],
                                                                        coordinateMesh[i + 1][j + 2],
                                                                        coordinateMesh[i][j + 2],
                                                                        coordinateMesh[i][j + 1],
                                                                        rightPart[element],
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

    double c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14, c15, c16, c17,
        c18, c19, c20, c21, c22, c23, c24, c25, c26, c27, c28, c29, c30, c31, c32, c33,
        c34, c35, c36, c37, c38, c39, c40, c41;

    setCoefficients(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14, c15, c16, c17,
                    c18, c19, c20, c21, c22, c23, c24, c25, c26, c27, c28, c29, c30, c31, c32, c33,
                    c34, c35, c36, c37, c38, c39, c40, c41,
                    A1, A2, A3, A4,
                    s1, s2, s3, s4, s5, s6, s7, s8,
                    z1, z2, z3, z4, z5);

    double *a = new double[8];
    double *b = new double[8];
    double *c = new double[8];
    double *d = new double[8];
    double *e = new double[8];
    double *f = new double[8];
    double *g = new double[8];
    double *t = new double[8];

    setFormFunctionsCoefficients(a, b, c, d, e, f, g, t, pointI, pointM);

    double valR1, valR2, valR3, valR4, valR5, valR6, valR7;

    for (int i = 0; i < 8; i++)
        for (int j = 0; j < 8; j++)
        {
            valR1 = c[j] * (c[i] * c1 + f[i] * c5 + d[i] * c10 + t[i] * c19 + g[i] * c24);
            valR2 = d[j] * (b[i] * c3 + c[i] * c11 + e[i] * c41 + f[i] * c13 + d[i] * c26 + t[i] * c30 + g[i] * c35);
            valR3 = f[j] * (c[i] * c6 + f[i] * c9 + d[i] * c14 + t[i] * c23 + g[i] * c33);
            valR4 = g[j] * (b[i] * c15 + c[i] * c25 + f[i] * c28 + e[i] * c31 + d[i] * c36 + t[i] * c38 + g[i] * c39);
            valR5 = t[j] * (b[i] * c7 + c[i] * c16 + e[i] * c20 + f[i] * c22 + d[i] * c29 + t[i] * c34 + g[i] * c37);
            valR6 = b[j] * (b[i] * c2 + d[i] * c4 + t[i] * c8 + e[i] * c12 + g[i] * c17);
            valR7 = e[j] * (b[i] * c40 + d[i] * c18 + t[i] * c21 + e[i] * c27 + g[i] * c32);
            localMatrix.setElement(i, j, valR1 + valR2 + valR3 + valR4 + valR5 + valR6 + valR7);
        }

    double zi = pointI.getX();
    double zk = pointK.getX();

    double si = pointI.getY();
    double sr = pointR.getY();

    double poisC1 = systemParameters.poissonC1;
    double poisC2 = systemParameters.poissonC2;

    for (int i = 0; i < 8; i++)
    {
        localRightPart.setElement(i,
                                  a[i] * poisC1 * poisC2 * z2 * (-sin(poisC2 * si) + sin(poisC2 * sr)) +
                                      b[i] * poisC1 * poisC2 * z3 * (-sin(poisC2 * si) + sin(poisC2 * sr)) +
                                      e[i] * poisC1 * poisC2 * z4 * (-sin(poisC2 * si) + sin(poisC2 * sr)) +
                                      c[i] * poisC1 * z2 * (-cos(poisC2 * si) + cos(poisC2 * sr) - poisC2 * si * sin(poisC2 * si) + poisC2 * sr * sin(poisC2 * sr)) +
                                      d[i] * poisC1 * z3 * (-cos(poisC2 * si) + cos(poisC2 * sr) - poisC2 * si * sin(poisC2 * si) + poisC2 * sr * sin(poisC2 * sr)) +
                                      g[i] * poisC1 * z4 * (-cos(poisC2 * si) + cos(poisC2 * sr) - poisC2 * si * sin(poisC2 * si) + poisC2 * sr * sin(poisC2 * sr)) + (1.0 / poisC2) * f[i] * poisC1 * z2 * (-2.0 * poisC2 * si * cos(poisC2 * si) + 2.0 * poisC2 * sr * cos(poisC2 * sr) + (2.0 - poisC2 * poisC2 * si * si) * sin(poisC2 * si) + (-2.0 + poisC2 * poisC2 * sr * sr) * sin(poisC2 * sr)) +
                                      (1.0 / poisC2) * poisC1 * t[i] * z3 * (-2.0 * poisC2 * si * cos(poisC2 * si) + 2.0 * poisC2 * sr * cos(poisC2 * sr) + (2.0 - poisC2 * poisC2 * si * si) * sin(poisC2 * si) + (-2.0 + poisC2 * poisC2 * sr * sr) * sin(poisC2 * sr)));
    }

    delete[] a;
    delete[] b;
    delete[] c;
    delete[] d;
    delete[] e;
    delete[] f;
    delete[] g;
    delete[] t;
}

void createGlobalPressureMatrixForRectangleElementSecondOrder(double **&matrixPressure,
                                                              SecondOrderRectangleContributionMatrix *&contributionMatrix,
                                                              double *&rightPartParam,
                                                              SecondOrderRectangleRightPart *&localRightPartsParam,
                                                              int n)
{
    int *globalNodeNumbersIJK = new int[8];
    int finiteElementNumber = 0;
    int numberOfFE = ((n - 1) / 2); //in one row

    int nSkippedRows = 0;
    int firstNodeInSmallRow;
    int nodeK;            //8, 16, ... second row in local element and beginning of the row
    int elementInRow = 0; //while iterating in row store here current FE

    for (int i = 0; i < n - 2; i += 2)
    {
        firstNodeInSmallRow = (i + 1) * n - nSkippedRows * numberOfFE;
        nodeK = (i + 2) * n - (nSkippedRows + 1) * numberOfFE;

        elementInRow = 0;
        for (int j = 0; j < n - 2; j += 2)
        {
            //for top triangle
            /*i*/ globalNodeNumbersIJK[0] = i * n + j - nSkippedRows * numberOfFE;
            /*j*/ globalNodeNumbersIJK[1] = firstNodeInSmallRow + elementInRow;
            /*k*/ globalNodeNumbersIJK[2] = nodeK + j;
            /*l*/ globalNodeNumbersIJK[3] = globalNodeNumbersIJK[2] + 1;
            /*m*/ globalNodeNumbersIJK[4] = globalNodeNumbersIJK[3] + 1;
            /*n*/ globalNodeNumbersIJK[5] = firstNodeInSmallRow + elementInRow + 1;
            /*r*/ globalNodeNumbersIJK[6] = globalNodeNumbersIJK[0] + 2;
            /*q*/ globalNodeNumbersIJK[7] = globalNodeNumbersIJK[0] + 1;

            for (int i1 = 0; i1 < 8; i1++)
                for (int i2 = 0; i2 < 8; i2++)
                    matrixPressure[globalNodeNumbersIJK[i1]][globalNodeNumbersIJK[i2]] +=
                        contributionMatrix[finiteElementNumber].getElement(i1, i2);

            for (int i1 = 0; i1 < 8; i1++)
                rightPartParam[globalNodeNumbersIJK[i1]] += localRightPartsParam[finiteElementNumber].getElement(i1);

            finiteElementNumber++;
            elementInRow++;
        }
        nSkippedRows++;
    }
}