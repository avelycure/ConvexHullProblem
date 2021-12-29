#include "SecondOrderTriangleFE.hpp"

void solveWithSecondOrderTriangleFE(TriangleContributionMatrixSecondOrder *&contributionMatrix,
                                    TriangleRightPartSecondOrder *&localRigthParts,
                                    Point **&coordinateMesh,
                                    double **&matrixPressure,
                                    double *&rightPart,
                                    SystemParameters &systemParameters)
{
    const int MATRIX_PRESSURE_SIZE = systemParameters.n * systemParameters.n;
    const int MATRIX_CONTRIBUTION_SIZE = (systemParameters.n - 1) * (systemParameters.n - 1) / 2;

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

    addBorderConditionsOnPressureValues(systemParameters,
                                        matrixPressure,
                                        rightPart,
                                        systemParameters.n,
                                        MATRIX_PRESSURE_SIZE);

    outputPressureMatrix(matrixPressure, MATRIX_PRESSURE_SIZE);
}

void createLocalMatrixForEveryElementQuadraticTriangles(TriangleContributionMatrixSecondOrder *&contributionMatrixParam,
                                                        Point **&coordinateMeshParam,
                                                        TriangleRightPartSecondOrder *&rightPartParam,
                                                        SystemParameters &systemParameters)
{
    int n = systemParameters.n;
    int finiteElementNumber = -1;

    for (int i = 0; i < n - 2; i += 2)
        for (int j = 0; j < n - 2; j += 2)
        {
            finiteElementNumber++;
            createLocalContributionMatrixForQuardaticTriangle(contributionMatrixParam[finiteElementNumber],
                                                              coordinateMeshParam[i][j],
                                                              coordinateMeshParam[i + 1][j],
                                                              coordinateMeshParam[i + 2][j],
                                                              coordinateMeshParam[i + 1][j + 1],
                                                              coordinateMeshParam[i][j + 2],
                                                              coordinateMeshParam[i][j + 1],
                                                              rightPartParam[finiteElementNumber],
                                                              systemParameters);
            finiteElementNumber++;
            createLocalContributionMatrixForQuardaticTriangle(contributionMatrixParam[finiteElementNumber],
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

void createLocalContributionMatrixForQuardaticTriangle(TriangleContributionMatrixSecondOrder localMatrix,
                                                       Point pointI,
                                                       Point pointJ,
                                                       Point pointK,
                                                       Point pointL,
                                                       Point pointM,
                                                       Point pointN,
                                                       TriangleRightPartSecondOrder localRightPart,
                                                       SystemParameters &systemParameters)
{
    double hMin = systemParameters.hMin;
    double k = systemParameters.k;

    double k1 = (pointM.getX() - pointK.getX()) / (pointM.getY() - pointK.getY());
    double k2 = pointK.getX() - pointK.getY() * ((pointM.getX() - pointK.getX()) / (pointM.getY() - pointK.getY()));

    double A1 = pow(hMin, 3.0);
    double A2 = 3.0 * pow(hMin, 2.0) * k;
    double A3 = 3.0 * hMin * pow(k, 2.0);
    double A4 = pow(k, 3.0);

    double s1 = (1.0 / 1.0) * (pointM.getY() - pointI.getY());
    double s2 = (1.0 / 2.0) * (pow(pointM.getY(), 2.0) - pow(pointI.getY(), 2.0));
    double s3 = (1.0 / 3.0) * (pow(pointM.getY(), 3.0) - pow(pointI.getY(), 3.0));
    double s4 = (1.0 / 4.0) * (pow(pointM.getY(), 4.0) - pow(pointI.getY(), 4.0));
    double s5 = (1.0 / 5.0) * (pow(pointM.getY(), 5.0) - pow(pointI.getY(), 5.0));
    double s6 = (1.0 / 6.0) * (pow(pointM.getY(), 6.0) - pow(pointI.getY(), 6.0));
    double s7 = (1.0 / 7.0) * (pow(pointM.getY(), 7.0) - pow(pointI.getY(), 7.0));
    double s8 = (1.0 / 8.0) * (pow(pointM.getY(), 8.0) - pow(pointI.getY(), 8.0));

    double zi = pointI.getX();

    double c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14, c15, c16, c18, c19;

    setCoefficients(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14, c15, c16, c18, c19,
                    s1, s2, s3, s4, s5, s6, s7,
                    k1, k2,
                    A1, A2, A3, A4,
                    zi);

    double R6 = 6.0 * systemParameters.mu * systemParameters.L * systemParameters.U * k /
                (systemParameters.Hn * systemParameters.Hn * systemParameters.pMin);

    double *a = new double[6];
    double *b = new double[6];
    double *c = new double[6];
    double *d = new double[6];
    double *e = new double[6];
    double *f = new double[6];

    setFormFunctionsCoefficients(a, b, c, d, e, f, pointI, pointJ, pointK, pointM, pointN);

    double valR1, valR2, valR3, valR4, valR5;
    for (int i = 0; i < 6; i++)
        for (int j = 0; j < 6; j++)
        {
            valR1 = c[j] * (c[i] * c1 + f[i] * c3 + e[i] * c13);
            valR2 = e[j] * (b[i] * c5 + d[i] * c8 + f[i] * c9 + c[i] * c14 + e[i] * c18);
            valR3 = f[j] * (c[i] * c4 + f[i] * c7 + e[i] * c11);
            valR4 = b[j] * (b[i] * c2 + e[i] * c6 + d[i] * c15);
            valR5 = d[j] * (e[i] * c10 + b[i] * c16 + d[i] * c19);
            localMatrix.setElement(i, j, valR1 + valR2 + valR3 + valR4 + valR5);
        }

    double poisC1 = systemParameters.poissonC1;
    double poisC2 = systemParameters.poissonC2;

    double rpc1, rpc2, rpc3, rpc4;

    double si = pointI.getY();
    double sr = pointM.getY();

    std::cout << "*******" << std::endl;
    std::cout << "ZERO: " << std::endl;
    std::cout << a[0] << " " << b[0] << " " << c[0] << " " << d[0] << " " << e[0] << " " << f[0] << std::endl;
    std::cout << "k1 + k2:" << k1 << " " << k2 << std::endl;
    std::cout << "si + zi:" << si << " " << zi << std::endl;
    std::cout << "sr + si:" << sr << " " << si << std::endl;
    std::cout << "*******" << std::endl;

    for (int i = 0; i < 6; i++)
    {
        rpc1 = -(1.0 / (6.0 * poisC2 * poisC2)) * poisC1 * cos(poisC2 * si) * (3.0 * c[i] * (4.0 * k1 * k2 * poisC2 * poisC2 * si + 3.0 * k1 * k1 * (-2.0 + poisC2 * poisC2 * si * si) + poisC2 * poisC2 * (k2 * k2 - zi * zi)) + 2.0 * (-18.0 * e[i] * k1 * k1 * k2 - 18.0 * d[i] * k1 * k1 * k1 * k2 + 3.0 * a[i] * k1 * k2 * poisC2 * poisC2 + e[i] * k2 * k2 * k2 * poisC2 * poisC2 + 3.0 * d[i] * k1 * k2 * k2 * k2 * poisC2 * poisC2 - 24.0 * e[i] * k1 * k1 * k1 * si - 18.0 * d[i] * k1 * k1 * k1 * k1 * si + 3.0 * a[i] * k1 * k1 * poisC2 * poisC2 * si + 6.0 * e[i] * k1 * k2 * k2 * poisC2 * poisC2 * si + 9.0 * d[i] * k1 * k1 * k2 * k2 * poisC2 * poisC2 * si + 9.0 * e[i] * k1 * k1 * k2 * poisC2 * poisC2 * si * si + 9.0 * d[i] * k1 * k1 * k1 * k2 * poisC2 * poisC2 * si * si + 4.0 * e[i] * k1 * k1 * k1 * poisC2 * poisC2 * si * si * si + 3.0 * d[i] * k1 * k1 * k1 * k1 * poisC2 * poisC2 * si * si * si + 3.0 * b[i] * k1 * (k2 * k2 * poisC2 * poisC2 + 2.0 * k1 * k2 * poisC2 * poisC2 * si + k1 * k1 * (-2.0 + poisC2 * poisC2 * si * si)) - e[i] * poisC2 * poisC2 * zi * zi * zi + 3.0 * f[i] * (2.0 * k1 * k1 * si * (-6.0 + poisC2 * poisC2 * si * si) + 3.0 * k1 * k2 * (-2.0 + poisC2 * poisC2 * si * si) + poisC2 * poisC2 * si * (k2 * k2 - zi * zi))));

        rpc2 = (1.0 / (6.0 * poisC2 * poisC2)) * poisC1 * cos(poisC2 * sr) * (3.0 * c[i] * (4.0 * k1 * k2 * poisC2 * poisC2 * sr + 3.0 * k1 * k1 * (-2.0 + poisC2 * poisC2 * sr * sr) + poisC2 * poisC2 * (k2 * k2 - zi * zi)) + 2.0 * (-18.0 * e[i] * k1 * k1 * k2 - 18.0 * d[i] * k1 * k1 * k1 * k2 + 3.0 * a[i] * k1 * k2 * poisC2 * poisC2 + e[i] * k2 * k2 * k2 * poisC2 * poisC2 + 3.0 * d[i] * k1 * k2 * k2 * k2 * poisC2 * poisC2 - 24.0 * e[i] * k1 * k1 * k1 * sr - 18.0 * d[i] * k1 * k1 * k1 * k1 * sr + 3.0 * a[i] * k1 * k1 * poisC2 * poisC2 * sr + 6.0 * e[i] * k1 * k2 * k2 * poisC2 * poisC2 * sr + 9.0 * d[i] * k1 * k1 * k2 * k2 * poisC2 * poisC2 * sr + 9.0 * e[i] * k1 * k1 * k2 * poisC2 * poisC2 * sr * sr + 9.0 * d[i] * k1 * k1 * k1 * k2 * poisC2 * poisC2 * sr * sr + 4.0 * e[i] * k1 * k1 * k1 * poisC2 * poisC2 * sr * sr * sr + 3.0 * d[i] * k1 * k1 * k1 * k1 * poisC2 * poisC2 * sr * sr * sr + 3.0 * b[i] * k1 * (k2 * k2 * poisC2 * poisC2 + 2.0 * k1 * k2 * poisC2 * poisC2 * sr + k1 * k1 * (-2.0 + poisC2 * poisC2 * sr * sr)) - e[i] * poisC2 * poisC2 * zi * zi * zi + 3.0 * f[i] * (2.0 * k1 * k1 * sr * (-6.0 + poisC2 * poisC2 * sr * sr) + 3.0 * k1 * k2 * (-2.0 + poisC2 * poisC2 * sr * sr) + poisC2 * poisC2 * sr * (k2 * k2 - zi * zi))));

        rpc3 = (-1.0 / (12.0 * poisC2 * poisC2 * poisC2)) * poisC1 * sin(poisC2 * si) * (72.0 * d[i] * k1 * k1 * k1 * k1 - 12.0 * a[i] * k1 * k1 * poisC2 * poisC2 - 24.0 * c[i] * k1 * k2 * poisC2 * poisC2 - 24.0 * b[i] * k1 * k1 * k2 * poisC2 * poisC2 - 36.0 * d[i] * k1 * k1 * k2 * k2 * poisC2 * poisC2 + 6.0 * a[i] * k2 * k2 * poisC2 * poisC2 * poisC2 * poisC2 + 4.0 * b[i] * k2 * k2 * k2 * poisC2 * poisC2 * poisC2 * poisC2 + 3.0 * d[i] * k2 * k2 * k2 * k2 * poisC2 * poisC2 * poisC2 * poisC2 - 36.0 * c[i] * k1 * k1 * poisC2 * poisC2 * si - 24.0 * b[i] * k1 * k1 * k1 * poisC2 * poisC2 * si - 72.0 * d[i] * k1 * k1 * k1 * k2 * poisC2 * poisC2 * si + 12.0 * a[i] * k1 * k2 * poisC2 * poisC2 * poisC2 * poisC2 * si + 6.0 * c[i] * k2 * k2 * poisC2 * poisC2 * poisC2 * poisC2 * si + 12.0 * b[i] * k1 * k2 * k2 * poisC2 * poisC2 * poisC2 * poisC2 * si + 12.0 * d[i] * k1 * k2 * k2 * k2 * poisC2 * poisC2 * poisC2 * poisC2 * si - 36.0 * d[i] * k1 * k1 * k1 * k1 * poisC2 * poisC2 * si * si + 6.0 * a[i] * k1 * k1 * poisC2 * poisC2 * poisC2 * poisC2 * si * si + 12.0 * c[i] * k1 * k2 * poisC2 * poisC2 * poisC2 * poisC2 * si * si + 12.0 * b[i] * k1 * k1 * k2 * poisC2 * poisC2 * poisC2 * poisC2 * si * si + 18.0 * d[i] * k1 * k1 * k2 * k2 * poisC2 * poisC2 * poisC2 * poisC2 * si * si + 6.0 * c[i] * k1 * k1 * poisC2 * poisC2 * poisC2 * poisC2 * si * si * si + 4.0 * b[i] * k1 * k1 * k1 * poisC2 * poisC2 * poisC2 * poisC2 * si * si * si + 12.0 * d[i] * k1 * k1 * k1 * k2 * poisC2 * poisC2 * poisC2 * poisC2 * si * si * si + 3.0 * d[i] * k1 * k1 * k1 * k1 * poisC2 * poisC2 * poisC2 * poisC2 * si * si * si * si -6.0 * a[i] * poisC2 * poisC2 * poisC2 * poisC2 * zi * zi - 6.0 * c[i] * poisC2 * poisC2 * poisC2 * poisC2 * si * zi * zi - 4.0 * b[i] * poisC2 * poisC2 * poisC2 * poisC2 * zi * zi * zi - 3.0 * d[i] * poisC2 * poisC2 * poisC2 * poisC2 * zi * zi * zi * zi + 6.0 * f[i] * (2.0 * k1 * k2 * poisC2 * poisC2 * si * (-6.0 + poisC2 * poisC2 * si * si) + k1 * k1 * (24.0 - 12.0 * poisC2 * poisC2 * si * si + poisC2 * poisC2 * poisC2 * poisC2 * si * si * si * si) + poisC2 * poisC2 * (-2.0 + poisC2 * poisC2 * si * si) * (k2 * k2 - zi * zi)) + 4.0 * e[i] * (3.0 * k1 * k1 * k2 * poisC2 * poisC2 * si * (-6.0 + poisC2 * poisC2 * si * si) + 3.0 * k1 * k2 * k2 * poisC2 * poisC2 * (-2 + poisC2 * poisC2 * si * si) + k1 * k1 * k1 * (24.0 - 12.0 * poisC2 * poisC2 * si * si + poisC2 * poisC2 * poisC2 * poisC2 * si * si * si * si) + poisC2 * poisC2 * poisC2 * poisC2 * si * (k2 * k2 * k2 - zi * zi * zi)));

        rpc4 = (1.0 / (12.0 * poisC2 * poisC2 * poisC2)) * poisC1 * sin(poisC2 * sr) * (72.0 * d[i] * k1 * k1 * k1 * k1 - 12.0 * a[i] * k1 * k1 * poisC2 * poisC2 - 24.0 * c[i] * k1 * k2 * poisC2 * poisC2 - 24.0 * b[i] * k1 * k1 * k2 * poisC2 * poisC2 - 36.0 * d[i] * k1 * k1 * k2 * k2 * poisC2 * poisC2 + 6.0 * a[i] * k2 * k2 * poisC2 * poisC2 * poisC2 * poisC2 + 4.0 * b[i] * k2 * k2 * k2 * poisC2 * poisC2 * poisC2 * poisC2 + 3.0 * d[i] * k2 * k2 * k2 * k2 * poisC2 * poisC2 * poisC2 * poisC2 - 36.0 * c[i] * k1 * k1 * poisC2 * poisC2 * sr - 24.0 * b[i] * k1 * k1 * k1 * poisC2 * poisC2 * sr - 72.0 * d[i] * k1 * k1 * k1 * k2 * poisC2 * poisC2 * sr + 12.0 * a[i] * k1 * k2 * poisC2 * poisC2 * poisC2 * poisC2 * sr + 6.0 * c[i] * k2 * k2 * poisC2 * poisC2 * poisC2 * poisC2 * sr + 12.0 * b[i] * k1 * k2 * k2 * poisC2 * poisC2 * poisC2 * poisC2 * sr + 12.0 * d[i] * k1 * k2 * k2 * k2 * poisC2 * poisC2 * poisC2 * poisC2 * sr - 36.0 * d[i] * k1 * k1 * k1 * k1 * poisC2 * poisC2 * sr * sr + 6.0 * a[i] * k1 * k1 * poisC2 * poisC2 * poisC2 * poisC2 * sr * sr + 12.0 * c[i] * k1 * k2 * poisC2 * poisC2 * poisC2 * poisC2 * sr * sr + 12.0 * b[i] * k1 * k1 * k2 * poisC2 * poisC2 * poisC2 * poisC2 * sr * sr + 18.0 * d[i] * k1 * k1 * k2 * k2 * poisC2 * poisC2 * poisC2 * poisC2 * sr * sr + 6.0 * c[i] * k1 * k1 * poisC2 * poisC2 * poisC2 * poisC2 * sr * sr * sr + 4.0 * b[i] * k1 * k1 * k1 * poisC2 * poisC2 * poisC2 * poisC2 * sr * sr * sr + 12.0 * d[i] * k1 * k1 * k1 * k2 * poisC2 * poisC2 * poisC2 * poisC2 * sr * sr * sr + 3.0 * d[i] * k1 * k1 * k1 * k1 * poisC2 * poisC2 * poisC2 * poisC2 * sr * sr * sr * sr - 6.0 * a[i] * poisC2 * poisC2 * poisC2 * poisC2 * zi * zi - 6.0 * c[i] * poisC2 * poisC2 * poisC2 * poisC2 * sr * zi * zi - 4.0 * b[i] * poisC2 * poisC2 * poisC2 * poisC2 * zi * zi * zi - 3.0 * d[i] * poisC2 * poisC2 * poisC2 * poisC2 * zi * zi * zi * zi + 6.0 * f[i] * (2.0 * k1 * k2 * poisC2 * poisC2 * sr * (-6.0 + poisC2 * poisC2 * sr * sr) + k1 * k1 * (24.0 - 12.0 * poisC2 * poisC2 * sr * sr + poisC2 * poisC2 * poisC2 * poisC2 * sr * sr * sr * sr) + poisC2 * poisC2 * (-2.0 + poisC2 * poisC2 * sr * sr) * (k2 * k2 - zi * zi)) + 4.0 * e[i] * (3.0 * k1 * k1 * k2 * poisC2 * poisC2 * sr * (-6.0 + poisC2 * poisC2 * sr * sr) + 3.0 * k1 * k2 * k2 * poisC2 * poisC2 * (-2.0 + poisC2 * poisC2 * sr * sr) + k1 * k1 * k1 * (24.0 - 12.0 * poisC2 * poisC2 * sr * sr + poisC2 * poisC2 * poisC2 * poisC2 * sr * sr * sr * sr) + poisC2 * poisC2 * poisC2 * poisC2 * sr * (k2 * k2 * k2 - zi * zi * zi)));

        if (i == 0)
            std::cout << "RPC: " << rpc1 << " " << rpc2 << " " << rpc3 << " "
                      << " " << rpc4 << std::endl;
        localRightPart.setElement(i, rpc1 + rpc2 + rpc3 + rpc4);
    }

    delete[] a;
    delete[] b;
    delete[] c;
    delete[] d;
    delete[] e;
    delete[] f;
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