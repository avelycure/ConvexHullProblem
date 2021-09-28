#include "header.hpp"
/*
*    LINEAR PART
*/
int solveWithHLinear(ContributionMatrix *&contributionMatrix,
                     RightPart *&localRigthParts,
                     Point **&coordinateMesh,
                     double **&matrixPressure,
                     double *&rightPart,
                     SystemPatemeters &systemParameters)
{
    const int MATRIX_PRESSURE_SIZE = systemParameters.n * systemParameters.n;
    const int MATRIX_CONTRIBUTION_SIZE = (systemParameters.n - 1) * (systemParameters.n - 1) * 2;

    initMatrix(matrixPressure, MATRIX_PRESSURE_SIZE, MATRIX_PRESSURE_SIZE);
    initMesh(coordinateMesh, systemParameters);

    initContributionMatrix(contributionMatrix, MATRIX_CONTRIBUTION_SIZE);
    initRightPart(localRigthParts, MATRIX_CONTRIBUTION_SIZE);
    initVector(rightPart, MATRIX_PRESSURE_SIZE);

    createLocalMatrixForEveryElementHLinear(contributionMatrix, coordinateMesh, localRigthParts, systemParameters);

    createGlobalPressureMatrixHLinear(matrixPressure, contributionMatrix, rightPart, localRigthParts, systemParameters.n);
    addBorderConditionsHLinear(matrixPressure, rightPart, systemParameters.n, MATRIX_PRESSURE_SIZE,
                               systemParameters.LOW_BORDER, systemParameters.HIGH_BORDER);

    outputPressureMatrix(matrixPressure, MATRIX_PRESSURE_SIZE);

    return 0;
}

int createLocalContributionMatrixForHLinearTop(ContributionMatrix localMatrix,
                                               Point pointI, Point pointJ, Point pointK,
                                               RightPart localRightPart,
                                               SystemPatemeters &systemParameters)
{
    double hMin = systemParameters.hMin;
    double k = systemParameters.k;

    // coefficients of line z = k1 * x + k2
    double k1 = (pointK.getX() - pointJ.getX()) / (pointK.getY() - pointJ.getY());
    double k2 = pointJ.getX() - pointJ.getY() * (pointK.getX() - pointJ.getX()) / (pointK.getY() - pointJ.getY());

    //coefficients of polynom
    double A1 = pow(hMin, 3.0);
    double A2 = 3.0 * pow(hMin, 2.0) * k;
    double A3 = 3.0 * hMin * pow(k, 2.0);
    double A4 = pow(k, 3.0);

    double s = pointK.getY() - pointI.getY();
    double s2 = 0.5 * (pow(pointK.getY(), 2.0) - pow(pointI.getY(), 2.0));
    double s3 = (1.0 / 3.0) * (pow(pointK.getY(), 3.0) - pow(pointI.getY(), 3.0));
    double s4 = 0.25 * (pow(pointK.getY(), 4.0) - pow(pointI.getY(), 4.0));
    double s5 = 0.2 * (pow(pointK.getY(), 5.0) - pow(pointI.getY(), 5.0));

    double zk = pointK.getX() - k2;

    //coefficients of T1
    double coefT1 = -k1 * A4 * s5 + (A4 * zk - A3 * k1) * s4 + (A3 * zk - A2 * k1) * s3 +
                    (-k1 * A1 + A2 * zk) * s2 + A1 * zk * s;

    //coefficients of T2
    double coefT2 = -k1 * s2 + zk * s;

    //coefficients of T3
    double coefT3 = -k1 * s3 + zk * s2;

    //coefficients of T4
    double coefT4 = -0.5 * k1 * k1 * s3 - k1 * k2 * s2 + 0.5 * (pointK.getX() * pointK.getX() - k2 * k2) * s;

    //coefficient R3
    double R3 = 6.0 * systemParameters.mu * systemParameters.L * systemParameters.U * k /
                (systemParameters.Hn * systemParameters.Hn * systemParameters.pMin);

    double valB, valG;
    double *a = new double[3];
    double *b = new double[3];
    double *c = new double[3];

    double area = countArea(pointI, pointJ, pointK);

    a[0] = pointJ.getX() * pointK.getY() - pointK.getX() * pointJ.getY();
    a[1] = pointK.getX() * pointI.getY() - pointI.getX() * pointK.getY();
    a[2] = pointI.getX() * pointJ.getY() - pointJ.getX() * pointI.getY();

    b[0] = pointJ.getY() - pointK.getY();
    b[1] = pointK.getY() - pointI.getY();
    b[2] = pointI.getY() - pointJ.getY();

    c[0] = pointK.getX() - pointJ.getX();
    c[1] = pointI.getX() - pointK.getX();
    c[2] = pointJ.getX() - pointI.getX();

    for (int i = 0; i < 3; i++)
    {
        a[i] = a[i] / (2.0 * area);
        b[i] = b[i] / (2.0 * area);
        c[i] = c[i] / (2.0 * area);
    }

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            valB = b[i] * b[j] * coefT1;
            valG = c[i] * c[j] * coefT1;
            localMatrix.setElement(i, j, valB + valG);
        }
    }

    //creating local vector
    for (int i = 0; i < 3; i++)
        localRightPart.vector[i] = -R3 * (coefT4 * c[i] + coefT3 * b[i] + coefT2 * a[i]);

    return 0;
}

int createLocalContributionMatrixForHLinearBottom(ContributionMatrix localMatrix,
                                                  Point pointI, Point pointJ, Point pointK,
                                                  RightPart localRightPart,
                                                  SystemPatemeters &systemParameters)
{
    // coefficients of line z = k1 * x + k2
    double k1 = (pointK.getX() - pointJ.getX()) / (pointK.getY() - pointJ.getY());
    double k2 = pointJ.getX() - pointJ.getY() * (pointK.getX() - pointJ.getX()) / (pointK.getY() - pointJ.getY());

    double hMin = systemParameters.hMin;
    double k = systemParameters.k;

    //coefficients of polynom
    double A1 = pow(hMin, 3.0);
    double A2 = 3.0 * pow(hMin, 2.0) * k;
    double A3 = 3.0 * hMin * pow(k, 2.0);
    double A4 = pow(k, 3.0);

    double s = pointI.getY() - pointK.getY();
    double s2 = 0.5 * (pow(pointI.getY(), 2.0) - pow(pointK.getY(), 2.0));
    double s3 = (1.0 / 3.0) * (pow(pointI.getY(), 3.0) - pow(pointK.getY(), 3.0));
    double s4 = 0.25 * (pow(pointI.getY(), 4.0) - pow(pointK.getY(), 4.0));
    double s5 = 0.2 * (pow(pointI.getY(), 5.0) - pow(pointK.getY(), 5.0));

    double zk = k2 - pointK.getX();

    //coefficients of T1
    double coefT1 = k1 * A4 * s5 + (A4 * zk + A3 * k1) * s4 + (A3 * zk + A2 * k1) * s3 +
                    (k1 * A1 + A2 * zk) * s2 + A1 * zk * s;

    //coefficients of T2
    double coefT2 = k1 * s2 + zk * s;

    //coefficients of T3
    double coefT3 = k1 * s3 + zk * s2;

    //coefficients of T4
    double coefT4 = 0.5 * k1 * k1 * s3 + k1 * k2 * s2 + 0.5 * (k2 * k2 - pointK.getX() * pointK.getX()) * s;

    //coefficient R3 6 * mu * k * L * ...
    double R3 = 6.0 * systemParameters.mu * systemParameters.L * systemParameters.U * k /
                (systemParameters.Hn * systemParameters.Hn * systemParameters.pMin);

    double valB, valG;
    double *a = new double[3];
    double *b = new double[3];
    double *c = new double[3];

    double area = countArea(pointI, pointJ, pointK);

    a[0] = pointJ.getX() * pointK.getY() - pointK.getX() * pointJ.getY();
    a[1] = pointK.getX() * pointI.getY() - pointI.getX() * pointK.getY();
    a[2] = pointI.getX() * pointJ.getY() - pointJ.getX() * pointI.getY();

    b[0] = pointJ.getY() - pointK.getY();
    b[1] = pointK.getY() - pointI.getY();
    b[2] = pointI.getY() - pointJ.getY();

    c[0] = pointK.getX() - pointJ.getX();
    c[1] = pointI.getX() - pointK.getX();
    c[2] = pointJ.getX() - pointI.getX();

    for (int i = 0; i < 3; i++)
    {
        a[i] = a[i] / (2.0 * area);
        b[i] = b[i] / (2.0 * area);
        c[i] = c[i] / (2.0 * area);
    }

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            valB = b[i] * b[j] * coefT1;
            valG = c[i] * c[j] * coefT1;
            localMatrix.setElement(i, j, valB + valG);
        }
    }

    //creating local vector
    for (int i = 0; i < 3; i++)
        localRightPart.vector[i] = -R3 * (coefT4 * c[i] + coefT3 * b[i] + coefT2 * a[i]);

    return 0;
}

void createLocalMatrixForEveryElementHLinear(ContributionMatrix *&contributionMatrixParam,
                                             Point **&coordinateMeshParam,
                                             RightPart *&rightPartParam,
                                             SystemPatemeters &systemParameters)
{
    double k = systemParameters.k;
    double hMin = systemParameters.hMin;
    int n = systemParameters.n;

    int finiteElementNumber = -1;
    for (int i = 0; i < n - 1; i++)
    {
        for (int j = 0; j < n - 1; j++)
        {
            finiteElementNumber++;
            createLocalContributionMatrixForHLinearTop(contributionMatrixParam[finiteElementNumber],
                                                       coordinateMeshParam[i][j],
                                                       coordinateMeshParam[i + 1][j],
                                                       coordinateMeshParam[i][j + 1],
                                                       rightPartParam[finiteElementNumber],
                                                       systemParameters);
            finiteElementNumber++;
            createLocalContributionMatrixForHLinearBottom(contributionMatrixParam[finiteElementNumber],
                                                          coordinateMeshParam[i + 1][j + 1],
                                                          coordinateMeshParam[i][j + 1],
                                                          coordinateMeshParam[i + 1][j],
                                                          rightPartParam[finiteElementNumber],
                                                          systemParameters);
        }
    }
}

void createGlobalPressureMatrixHLinear(double **&matrixPressure, ContributionMatrix *&contributionMatrix,
                                       double *&rightPartParam, RightPart *&localRightPartsParam, int n)
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
            {
                for (int iterator2 = 0; iterator2 < 3; iterator2++)
                {
                    matrixPressure[globalNodeNumbersIJK[iterator1]][globalNodeNumbersIJK[iterator2]] +=
                        contributionMatrix[finiteElementNumber].matrix[iterator1][iterator2];
                }
            }

            for (int iterator1 = 0; iterator1 < 3; iterator1++)
            {
                rightPartParam[globalNodeNumbersIJK[iterator1]] += localRightPartsParam[finiteElementNumber].vector[iterator1];
            }
            finiteElementNumber++;

            //for bottom triangle
            globalNodeNumbersIJK[0] = globalNodeNumbersIJK[1] + 1.0;
            swap(globalNodeNumbersIJK[1], globalNodeNumbersIJK[2]);

            for (int iterator1 = 0; iterator1 < 3; iterator1++)
            {
                for (int iterator2 = 0; iterator2 < 3; iterator2++)
                {
                    matrixPressure[globalNodeNumbersIJK[iterator1]][globalNodeNumbersIJK[iterator2]] +=
                        contributionMatrix[finiteElementNumber].matrix[iterator1][iterator2];
                }
            }

            for (int iterator1 = 0; iterator1 < 3; iterator1++)
            {
                rightPartParam[globalNodeNumbersIJK[iterator1]] += localRightPartsParam[finiteElementNumber].vector[iterator1];
            }
            finiteElementNumber++;
        }
    }
}

void addBorderConditionsHLinear(double **&matrixResult,
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

    fstream myFile;
    myFile.open("data/rightPart.txt", fstream::out);
    for (int i = 0; i < MATRIX_PRESSURE_SIZE; i++)
        myFile << rightPartParam[i] << endl;
}
