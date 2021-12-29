#include "FirstOrderRectangleFE.hpp"

void solveWithFirstOrderRectangleFE(FirstOrderRectangleContributionMatrix *&contributionMatrix,
                                    FirstOrderRectangleRightPart *&localRigthParts,
                                    Point **&coordinateMesh,
                                    double **&matrixPressure,
                                    double *&rightPart,
                                    SystemParameters &systemParameters)
{
    const int MATRIX_PRESSURE_SIZE = systemParameters.n * systemParameters.n;
    const int MATRIX_CONTRIBUTION_SIZE = (systemParameters.n - 1) * (systemParameters.n - 1);

    initMatrix(matrixPressure, MATRIX_PRESSURE_SIZE, MATRIX_PRESSURE_SIZE);
    for (int i = 0; i < MATRIX_PRESSURE_SIZE; i++)
        for (int j = 0; j < MATRIX_PRESSURE_SIZE; j++)
            matrixPressure[i][j] = 0.0;
    initMesh(coordinateMesh, systemParameters);

    contributionMatrix = new FirstOrderRectangleContributionMatrix[MATRIX_CONTRIBUTION_SIZE];

    localRigthParts = new FirstOrderRectangleRightPart[MATRIX_CONTRIBUTION_SIZE];

    rightPart = new double[MATRIX_PRESSURE_SIZE];
    for (int i = 0; i < MATRIX_PRESSURE_SIZE; i++)
        rightPart[i] = 0.0;

    createLocalMatrixForEveryRectangleElement(contributionMatrix,
                                              coordinateMesh,
                                              localRigthParts,
                                              systemParameters);

    createGlobalPressureMatrixForRectangleElement(matrixPressure,
                                                  contributionMatrix,
                                                  rightPart,
                                                  localRigthParts,
                                                  systemParameters.n);

    addBorderConditionsOnPressureValues(matrixPressure,
                                        rightPart,
                                        systemParameters.n,
                                        MATRIX_PRESSURE_SIZE,
                                        systemParameters.LOW_BORDER,
                                        systemParameters.HIGH_BORDER);

    outputPressureMatrix(matrixPressure, MATRIX_PRESSURE_SIZE);
}

void createLocalMatrixForEveryRectangleElement(FirstOrderRectangleContributionMatrix *&contributionMatrixParam,
                                               Point **&coordinateMeshParam,
                                               FirstOrderRectangleRightPart *&rightPartParam,
                                               SystemParameters &systemParameters)
{
    int n = systemParameters.n;
    int element = 0;

    for (int i = 0; i < n - 1; i++)
        for (int j = 0; j < n - 1; j++)
        {
            createLocalContributionMatrixForRectangleElement(contributionMatrixParam[element],
                                                             coordinateMeshParam[i][j],
                                                             coordinateMeshParam[i + 1][j],
                                                             coordinateMeshParam[i + 1][j + 1],
                                                             coordinateMeshParam[i][j + 1],
                                                             rightPartParam[element],
                                                             systemParameters);
            element++;
        }
}

void createLocalContributionMatrixForRectangleElement(FirstOrderRectangleContributionMatrix &localMatrix,
                                                      Point pointI,
                                                      Point pointJ,
                                                      Point pointK,
                                                      Point pointM,
                                                      FirstOrderRectangleRightPart &localRightPart,
                                                      SystemParameters &systemParameters)
{
    double hMin = systemParameters.hMin;
    double k = systemParameters.k;

    //coefficients of polynom
    double A1 = pow(hMin, 3.0);
    double A2 = 3.0 * pow(hMin, 2.0) * k;
    double A3 = 3.0 * hMin * pow(k, 2.0);
    double A4 = pow(k, 3.0);

    double s1 = pointM.getY() - pointI.getY();
    double s2 = (1.0 / 2.0) * (pow(pointM.getY(), 2.0) - pow(pointI.getY(), 2.0));
    double s3 = (1.0 / 3.0) * (pow(pointM.getY(), 3.0) - pow(pointI.getY(), 3.0));
    double s4 = (1.0 / 4.0) * (pow(pointM.getY(), 4.0) - pow(pointI.getY(), 4.0));
    double s5 = (1.0 / 5.0) * (pow(pointM.getY(), 5.0) - pow(pointI.getY(), 5.0));
    double s6 = (1.0 / 6.0) * (pow(pointM.getY(), 6.0) - pow(pointI.getY(), 6.0));

    double z1 = pointJ.getX() - pointI.getX();
    double z2 = (1.0 / 2.0) * (pow(pointJ.getX(), 2.0) - pow(pointI.getX(), 2.0));
    double z3 = (1.0 / 3.0) * (pow(pointJ.getX(), 3.0) - pow(pointI.getX(), 3.0));

    //before cl R1
    double c1 = A1 * s1 * z1 + A2 * s2 * z1 + A3 * s3 * z1 + A4 * s4 * z1;

    //before bl R2
    double c2 = A1 * s1 * z1 + A2 * s2 * z1 + A3 * s3 * z1 + A4 * s4 * z1;

    //before dl R2
    double c3 = A1 * s2 * z1 + A2 * s3 * z1 + A3 * s4 * z1 + A4 * s5 * z1;

    //before blR4
    double c4 = A1 * s2 * z1 + A2 * s3 * z1 + A3 * s4 * z1 + A4 * s5 * z1;

    //before dlR1
    double c5 = A1 * s1 * z2 + A2 * s2 * z2 + A3 * s3 * z2 + A4 * s4 * z2;

    //before clR4
    double c6 = A1 * s1 * z2 + A2 * s2 * z2 + A3 * s3 * z2 + A4 * s4 * z2;

    //before dlR4
    double c7 = A1 * s3 * z1 + A2 * s4 * z1 + A3 * s5 * z1 + A4 * s6 * z1 +
                A1 * s1 * z3 + A2 * s2 * z3 + A3 * s3 * z3 + A4 * s4 * z3;

    //coefficient R3 6 * mu * k * L * ...
    double R3 = 6.0 * systemParameters.mu * systemParameters.L * systemParameters.U * k /
                (systemParameters.Hn * systemParameters.Hn * systemParameters.pMin);

    double *a = new double[4];
    double *b = new double[4];
    double *c = new double[4];
    double *d = new double[4];

    double xi = pointI.getX();
    double yi = pointI.getY();

    double xk = pointK.getX();
    double yk = pointK.getY();

    a[0] = (xk * yk) / ((xi - xk) * (yi - yk));
    a[1] = (xi * yk) / ((xi - xk) * (-yi + yk));
    a[2] = (xi * yi) / ((xi - xk) * (yi - yk));
    a[3] = (xk * yi) / ((-xi + xk) * (yi - yk));

    b[0] = yk / ((xi - xk) * (-yi + yk));
    b[1] = yk / ((xi - xk) * (yi - yk));
    b[2] = yi / ((xi - xk) * (-yi + yk));
    b[3] = yi / ((xi - xk) * (yi - yk));

    c[0] = xk / ((-xi + xk) * (yi - yk));
    c[1] = xi / ((xi - xk) * (yi - yk));
    c[2] = xi / ((xi - xk) * (-yi + yk));
    c[3] = xk / ((xi - xk) * (yi - yk));

    d[0] = 1.0 / ((xi - xk) * (yi - yk));
    d[1] = 1.0 / ((-xi + xk) * (yi - yk));
    d[2] = 1.0 / ((xi - xk) * (yi - yk));
    d[3] = 1.0 / ((-xi + xk) * (yi - yk));

    double valR1, valR2, valR4;

    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
        {
            valR1 = c[j] * (c[i] * c1 + d[i] * c5);
            valR2 = b[j] * (b[i] * c2 + d[i] * c3);
            valR4 = d[j] * (b[i] * c4 + c[i] * c6 + d[i] * c7);
            localMatrix.setElement(i, j, valR1 + valR2 + valR4);
        }

    //creating local vector
    for (int i = 0; i < 4; i++)
        localRightPart.setElement(i, -R3 * (a[i] * s1 * z1 + c[i] * s2 * z1 + b[i] * s1 * z2 + d[i] * s2 * z2));

    delete[] a;
    delete[] b;
    delete[] c;
    delete[] d;
}

void createGlobalPressureMatrixForRectangleElement(double **&matrixPressure,
                                                   FirstOrderRectangleContributionMatrix *&contributionMatrix,
                                                   double *&rightPartParam,
                                                   FirstOrderRectangleRightPart *&localRightPartsParam,
                                                   int n)
{
    int *globalNodeNumbersIJK = new int[4];
    int finiteElementNumber = 0;

    for (int i = 0; i < n - 1; i++)
        for (int j = 0; j < n - 1; j++)
        {
            //for top triangle
            /*i*/ globalNodeNumbersIJK[0] = i * n + j;
            /*j*/ globalNodeNumbersIJK[1] = globalNodeNumbersIJK[0] + n;
            /*k*/ globalNodeNumbersIJK[2] = globalNodeNumbersIJK[0] + n + 1;
            /*m*/ globalNodeNumbersIJK[3] = globalNodeNumbersIJK[0] + 1;

            for (int i1 = 0; i1 < 4; i1++)
                for (int i2 = 0; i2 < 4; i2++)
                    matrixPressure[globalNodeNumbersIJK[i1]][globalNodeNumbersIJK[i2]] +=
                        contributionMatrix[finiteElementNumber].getElement(i1, i2);

            for (int i1 = 0; i1 < 4; i1++)
                rightPartParam[globalNodeNumbersIJK[i1]] += localRightPartsParam[finiteElementNumber].getElement(i1);

            finiteElementNumber++;
        }

    delete[] globalNodeNumbersIJK;
}