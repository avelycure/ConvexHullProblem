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
                                               systemParameters.LOW_BORDER,
                                               systemParameters.HIGH_BORDER);

    /*for (int i = 0; i < MATRIX_CONTRIBUTION_SIZE; i++)
    {
        displayMatrix(contributionMatrix[i].getMatrix(), 8, 8);
        std::cout << std::endl;
    }*/

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

    std::cout << "S: " << std::endl;
    std::cout << s1 << " " << s2 << " " << s3 << " " << s4 << " " << s5 << " " << s6 << " " << s7 << " " << s8 << std::endl;
    std::cout << "Z: " << std::endl;
    std::cout << z1 << " " << z2 << " " << z3 << " " << z4 << " " << z5 << std::endl;

    //clR1
    double c1 = A1 * s1 * z1 + A2 * s2 * z1 + A3 * s3 * z1 + A4 * s4 * z1;

    //blR6
    double c2 = A1 * s1 * z1 + A2 * s2 * z1 + A3 * s3 * z1 + A4 * s4 * z1;

    //blR2
    double c3 = A1 * s2 * z1 + A2 * s3 * z1 + A3 * s4 * z1 + A4 * s5 * z1;

    //dlR6
    double c4 = A1 * s2 * z1 + A2 * s3 * z1 + A3 * s4 * z1 + A4 * s5 * z1;

    //flR1
    double c5 = 2.0 * A1 * s2 * z1 + 2.0 * A2 * s3 * z1 + 2.0 * A3 * s4 * z1 + 2.0 * A4 * s5 * z1;

    //clR3
    double c6 = 2.0 * A1 * s2 * z1 + 2.0 * A2 * s3 * z1 + 2.0 * A3 * s4 * z1 + 2.0 * A4 * s5 * z1;

    //blR5
    double c7 = A1 * s3 * z1 + A2 * s4 * z1 + A3 * s5 * z1 + A4 * s6 * z1;

    //tlR6
    double c8 = A1 * s3 * z1 + A2 * s4 * z1 + A3 * s5 * z1 + A4 * s6 * z1;

    //flR3
    double c9 = 4.0 * A1 * s3 * z1 + 4.0 * A2 * s4 * z1 + 4.0 * A3 * s5 * z1 + 4.0 * A4 * s6 * z1;

    //dlR1
    double c10 = A1 * s1 * z2 + A2 * s2 * z2 + A3 * s3 * z2 + A4 * s4 * z2;

    //clR2
    double c11 = A1 * s1 * z2 + A2 * s2 * z2 + A3 * s3 * z2 + A4 * s4 * z2;

    //elR6
    double c12 = 2.0 * A1 * s1 * z2 + 2.0 * A2 * s2 * z2 + 2.0 * A3 * s3 * z2 + 2.0 * A4 * s4 * z2;

    //blR7
    double c40 = 2.0 * A1 * s1 * z2 + 2.0 * A2 * s2 * z2 + 2.0 * A3 * s3 * z2 + 2.0 * A4 * s4 * z2;

    //elR2
    double c41 = 2.0 * A1 * s2 * z2 + 2.0 * A2 * s3 * z2 + 2.0 * A3 * s4 * z2 + 2.0 * A4 * s5 * z2;

    //flR2
    double c13 = 2.0 * A1 * s2 * z2 + 2.0 * A2 * s3 * z2 + 2.0 * A3 * s4 * z2 + 2.0 * A4 * s5 * z2;

    //dlR3
    double c14 = 2.0 * A1 * s2 * z2 + 2.0 * A2 * s3 * z2 + 2.0 * A3 * s4 * z2 + 2.0 * A4 * s5 * z2;

    //blR4
    double c15 = 2.0 * A1 * s2 * z2 + 2.0 * A2 * s3 * z2 + 2.0 * A3 * s4 * z2 + 2.0 * A4 * s5 * z2;

    //clR5
    double c16 = 2.0 * A1 * s2 * z2 + 2.0 * A2 * s3 * z2 + 2.0 * A3 * s4 * z2 + 2.0 * A4 * s5 * z2;

    //glR6
    double c17 = 2.0 * A1 * s2 * z2 + 2.0 * A2 * s3 * z2 + 2.0 * A3 * s4 * z2 + 2.0 * A4 * s5 * z2;

    //dlR7
    double c18 = 2 * A1 * s2 * z2 + 2 * A2 * s3 * z2 + 2 * A3 * s4 * z2 + 2 * A4 * s5 * z2;

    //tlR1
    double c19 = 2.0 * A1 * s2 * z2 + 2.0 * A2 * s3 * z2 + 2.0 * A3 * s4 * z2 + 2.0 * A4 * s5 * z2;

    //elR5
    double c20 = 2.0 * A1 * s3 * z2 + 2.0 * A2 * s4 * z2 + 2.0 * A3 * s5 * z2 + 2.0 * A4 * s6 * z2;

    //tlR7
    double c21 = 2.0 * A1 * s3 * z2 + 2.0 * A2 * s4 * z2 + 2.0 * A3 * s5 * z2 + 2.0 * A4 * s6 * z2;

    //flR5
    double c22 = 4.0 * A1 * s3 * z2 + 4.0 * A2 * s4 * z2 + 4.0 * A3 * s5 * z2 + 4.0 * A4 * s6 * z2;

    //tlR3
    double c23 = 4.0 * A1 * s3 * z2 + 4.0 * A2 * s4 * z2 + 4.0 * A3 * s5 * z2 + 4.0 * A4 * s6 * z2;

    //glR1
    double c24 = A1 * s1 * z3 + A2 * s2 * z3 + A3 * s3 * z3 + A4 * s4 * z3;

    //clR4
    double c25 = A1 * s1 * z3 + A2 * s2 * z3 + A3 * s3 * z3 + A4 * s4 * z3;

    //dlR2
    double c26 = A1 * s3 * z1 + A2 * s4 * z1 + A3 * s5 * z1 + A4 * s6 * z1 +
                 A1 * s1 * z3 + A2 * s2 * z3 + A3 * s3 * z3 + A4 * s4 * z3;

    //elR7
    double c27 = 4.0 * A1 * s1 * z3 + 4.0 * A2 * s2 * z3 + 4.0 * A3 * s3 * z3 + 4.0 * A4 * s4 * z3;

    //flR4
    double c28 = 2.0 * A1 * s2 * z3 + 2.0 * A2 * s3 * z3 + 2.0 * A3 * s4 * z3 + 2.0 * A4 * s5 * z3;

    //dlR5
    double c29 = A1 * s4 * z1 + A2 * s5 * z1 + A3 * s6 * z1 + A4 * s7 * z1 +
                 2.0 * A1 * s2 * z3 + 2.0 * A2 * s3 * z3 + 2.0 * A3 * s4 * z3 + 2.0 * A4 * s5 * z3;

    //tlR2
    double c30 = A1 * s4 * z1 + A2 * s5 * z1 + A3 * s6 * z1 + A4 * s7 * z1 +
                 2.0 * A1 * s2 * z3 + 2.0 * A2 * s3 * z3 + 2.0 * A3 * s4 * z3 + 2.0 * A4 * s5 * z3;

    //elR4
    double c31 = 4.0 * A1 * s2 * z3 + 4.0 * A2 * s3 * z3 + 4.0 * A3 * s4 * z3 + 4.0 * A4 * s5 * z3;

    //glR7
    double c32 = 4.0 * A1 * s2 * z3 + 4.0 * A2 * s3 * z3 + 4.0 * A3 * s4 * z3 + 4.0 * A4 * s5 * z3;

    //glR3
    double c33 = 2.0 * A1 * s2 * z3 + 2.0 * A2 * s3 * z3 + 2.0 * A3 * s4 * z3 +
                 2.0 * A4 * s5 * z3;

    //tlR5
    double c34 = A1 * s5 * z1 + A2 * s6 * z1 + A3 * s7 * z1 + A4 * s8 * z1 +
                 4.0 * A1 * s3 * z3 + 4.0 * A2 * s4 * z3 + 4.0 * A3 * s5 * z3 + 4.0 * A4 * s6 * z3;

    //glR2
    double c35 = 2.0 * A1 * s3 * z2 + 2.0 * A2 * s4 * z2 + 2.0 * A3 * s5 * z2 + 2.0 * A4 * s6 * z2 +
                 A1 * s1 * z4 + A2 * s2 * z4 + A3 * s3 * z4 + A4 * s4 * z4;

    //dlR4
    double c36 = 2.0 * A1 * s3 * z2 + 2.0 * A2 * s4 * z2 + 2.0 * A3 * s5 * z2 + 2.0 * A4 * s6 * z2 +
                 A1 * s1 * z4 + A2 * s2 * z4 + A3 * s3 * z4 + A4 * s4 * z4;

    //glR5
    double c37 = 2.0 * A1 * s4 * z2 + 2.0 * A2 * s5 * z2 + 2.0 * A3 * s6 * z2 + 2.0 * A4 * s7 * z2 +
                 2.0 * A1 * s2 * z4 + 2.0 * A2 * s3 * z4 + 2.0 * A3 * s4 * z4 + 2.0 * A4 * s5 * z4;

    //tlR4
    double c38 = 2.0 * A1 * s4 * z2 + 2.0 * A2 * s5 * z2 + 2.0 * A3 * s6 * z2 + 2.0 * A4 * s7 * z2 +
                 2.0 * A1 * s2 * z4 + 2.0 * A2 * s3 * z4 + 2.0 * A3 * s4 * z4 + 2.0 * A4 * s5 * z4;

    //glR4
    double c39 = 4.0 * A1 * s3 * z3 + 4.0 * A2 * s4 * z3 + 4.0 * A3 * s5 * z3 + 4.0 * A4 * s6 * z3 +
                 A1 * s1 * z5 + A2 * s2 * z5 + A3 * s3 * z5 + A4 * s4 * z5;

    //coefficient R8 = 6 * mu * k * L * ...
    double R8 = 6.0 * systemParameters.mu * systemParameters.L * systemParameters.U * k /
                (systemParameters.Hn * systemParameters.Hn * systemParameters.pMin);

    double *a = new double[8];
    double *b = new double[8];
    double *c = new double[8];
    double *d = new double[8];
    double *e = new double[8];
    double *f = new double[8];
    double *g = new double[8];
    double *t = new double[8];

    double xm = pointM.getX();
    double ym = pointM.getY();

    double xi = pointI.getX();
    double yi = pointI.getY();

    a[0] = (xm * ym * (xi * (-3 * yi + ym) + xm * (yi + ym))) / ((xi - xm) * (xi - xm) * (yi - ym) * (yi - ym));
    a[1] = (4 * xi * xm * ym) / ((xi - xm) * (xi - xm) * (yi - ym));
    a[2] = (xi * ym * (xm * (-3 * yi + ym) + xi * (yi + ym))) / ((xi - xm) * (xi - xm) * (yi - ym) * (yi - ym));
    a[3] = -((4 * xi * yi * ym) / ((xi - xm) * (yi - ym) * (yi - ym)));
    a[4] = (xi * yi * (xm * (yi - 3 * ym) + xi * (yi + ym))) / ((xi - xm) * (xi - xm) * (yi - ym) * (yi - ym));
    a[5] = (4 * xi * xm * yi) / ((xi - xm) * (xi - xm) * (-yi + ym));
    a[6] = (xm * yi * (xi * (yi - 3 * ym) + xm * (yi + ym))) / ((xi - xm) * (xi - xm) * (yi - ym) * (yi - ym));
    a[7] = (4 * xm * yi * ym) / ((xi - xm) * (yi - ym) * (yi - ym));

    b[0] = (ym * ((3 * xi + xm) * yi - (xi + 3 * xm) * ym)) / ((xi - xm) * (xi - xm) * (yi - ym) * (yi - ym));
    b[1] = -((4 * (xi + xm) * ym) / ((xi - xm) * (xi - xm) * (yi - ym)));
    b[2] = (ym * ((xi + 3 * xm) * yi - (3 * xi + xm) * ym)) / ((xi - xm) * (xi - xm) * (yi - ym) * (yi - ym));
    b[3] = (4 * yi * ym) / ((xi - xm) * (yi - ym) * (yi - ym));
    b[4] = (yi * (-(3 * xi + xm) * yi + (xi + 3 * xm) * ym)) / ((xi - xm) * (xi - xm) * (yi - ym) * (yi - ym));
    b[5] = (4 * (xi + xm) * yi) / ((xi - xm) * (xi - xm) * (yi - ym));
    b[6] = (yi * (-(xi + 3 * xm) * yi + (3 * xi + xm) * ym)) / ((xi - xm) * (xi - xm) * (yi - ym) * (yi - ym));
    b[7] = -((4 * yi * ym) / ((xi - xm) * (yi - ym) * (yi - ym)));

    c[0] = (xm * (xi * (3 * yi + ym) - xm * (yi + 3 * ym))) / ((xi - xm) * (xi - xm) * (yi - ym) * (yi - ym));
    c[1] = (4 * xi * xm) / ((xi - xm) * (xi - xm) * (-yi + ym));
    c[2] = (xi * (xm * (3 * yi + ym) - xi * (yi + 3 * ym))) / ((xi - xm) * (xi - xm) * (yi - ym) * (yi - ym));
    c[3] = (4 * xi * (yi + ym)) / ((xi - xm) * (yi - ym) * (yi - ym));
    c[4] = (xi * (-xi * (3 * yi + ym) + xm * (yi + 3 * ym))) / ((xi - xm) * (xi - xm) * (yi - ym) * (yi - ym));
    c[5] = (4 * xi * xm) / ((xi - xm) * (xi - xm) * (yi - ym));
    c[6] = (xm * (-xm * (3 * yi + ym) + xi * (yi + 3 * ym))) / ((xi - xm) * (xi - xm) * (yi - ym) * (yi - ym));
    c[7] = -((4 * xm * (yi + ym)) / ((xi - xm) * (yi - ym) * (yi - ym)));

    d[0] = -((xm * (yi - 5 * ym) + xi * (3 * yi + ym)) / ((xi - xm) * (xi - xm) * (yi - ym) * (yi - ym)));
    d[1] = (4 * (xi + xm)) / ((xi - xm) * (xi - xm) * (yi - ym));
    d[2] = -((xi * (yi - 5 * ym) + xm * (3 * yi + ym)) / ((xi - xm) * (xi - xm) * (yi - ym) * (yi - ym)));
    d[3] = -((4 * (yi + ym)) / ((xi - xm) * (yi - ym) * (yi - ym)));
    d[4] = -((xi * (-5 * yi + ym) + xm * (yi + 3 * ym)) / ((xi - xm) * (xi - xm) * (yi - ym) * (yi - ym)));
    d[5] = (4 * (xi + xm)) / ((xi - xm) * (xi - xm) * (-yi + ym));
    d[6] = -((xm * (-5 * yi + ym) + xi * (yi + 3 * ym)) / ((xi - xm) * (xi - xm) * (yi - ym) * (yi - ym)));
    d[7] = (4 * (yi + ym)) / ((xi - xm) * (yi - ym) * (yi - ym));

    e[0] = (2 * ym) / ((xi - xm) * (xi - xm) * (-yi + ym));
    e[1] = (4 * ym) / ((xi - xm) * (xi - xm) * (yi - ym));
    e[2] = (2 * ym) / ((xi - xm) * (xi - xm) * (-yi + ym));
    e[3] = 0.0;
    e[4] = (2 * yi) / ((xi - xm) * (xi - xm) * (yi - ym));
    e[5] = (4 * yi) / ((xi - xm) * (xi - xm) * (-yi + ym));
    e[6] = (2 * yi) / ((xi - xm) * (xi - xm) * (yi - ym));
    e[7] = 0;

    f[0] = (2 * xm) / ((-xi + xm) * (yi - ym) * (yi - ym));
    f[1] = 0.0;
    f[2] = (2 * xi) / ((xi - xm) * (yi - ym) * (yi - ym));
    f[3] = -((4 * xi) / ((xi - xm) * (yi - ym) * (yi - ym)));
    f[4] = (2 * xi) / ((xi - xm) * (yi - ym) * (yi - ym));
    f[5] = 0.0;
    f[6] = (2 * xm) / ((-xi + xm) * (yi - ym) * (yi - ym));
    f[7] = (4 * xm) / ((xi - xm) * (yi - ym) * (yi - ym));

    g[0] = 2 / ((xi - xm) * (xi - xm) * (yi - ym));
    g[1] = 4 / ((xi - xm) * (xi - xm) * (-yi + ym));
    g[2] = 2 / ((xi - xm) * (xi - xm) * (yi - ym));
    g[3] = 0.0;
    g[4] = 2 / ((xi - xm) * (xi - xm) * (-yi + ym));
    g[5] = 4 / ((xi - xm) * (xi - xm) * (yi - ym));
    g[6] = 2 / ((xi - xm) * (xi - xm) * (-yi + ym));
    g[7] = 0.0;

    t[0] = 2 / ((xi - xm) * (yi - ym) * (yi - ym));
    t[1] = 0.0;
    t[2] = -(2 / ((xi - xm) * (yi - ym) * (yi - ym)));
    t[3] = 4 / ((xi - xm) * (yi - ym) * (yi - ym));
    t[4] = -(2 / ((xi - xm) * (yi - ym) * (yi - ym)));
    t[5] = 0.0;
    t[6] = 2 / ((xi - xm) * (yi - ym) * (yi - ym));
    t[7] = -(4 / ((xi - xm) * (yi - ym) * (yi - ym)));
    //std::cout << "Coefs: (" << pointI.getX() << ", " <<pointI.getY() << ") ("<< pointM.getX() << ", " <<pointM.getY() << ") "<<std::endl;

    //for (int i = 0; i < 8; i++)
    //    std::cout << a[i] << " " << b[i] << " " << c[i] << " " << d[i] << " " << e[i] << " " << f[i] << " " << g[i] << " " << t[i] << std::endl;
    std::cout << "********" << std::endl;
    std::cout << c1 << " " << c2 << " " << c3 << " " << c4 << " " << c5 << " " << c6 << " " << c7 << " " << c8 << std::endl;
    std::cout << c9 << " " << c10 << " " << c11 << " " << c12 << " " << c13 << " " << c14 << " " << c15 << " " << c16 << std::endl;
    std::cout << c17 << " " << c18 << " " << c19 << " " << c20 << " " << c21 << " " << c22 << " " << c23 << " " << c24 << std::endl;
    std::cout << c25 << " " << c26 << " " << c27 << " " << c28 << " " << c29 << " " << c30 << " " << c31 << " " << c32 << std::endl;
    std::cout << c33 << " " << c34 << " " << c35 << " " << c36 << " " << c37 << " " << c38 << " " << c39 << " " << c40 << std::endl;
    std::cout << c41 << std::endl;
    std::cout << "********" << std::endl;

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

    for (int i = 0; i < 8; i++)
        localRightPart.setElement(i, -R8 * (a[i] * s1 * z1 + c[i] * s2 * z1 + f[i] * s3 * z1 + b[i] * s1 * z2 +
                                            d[i] * s2 * z2 + t[i] * s3 * z2 + e[i] * s1 * z3 + g[i] * s2 * z3));

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