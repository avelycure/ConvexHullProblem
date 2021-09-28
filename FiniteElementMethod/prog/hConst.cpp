#include "header.hpp"

/**
 * Constant height. Special case of linear changing height
 * */
void solveWithHConst(ContributionMatrix *&contributionMatrix,
                     Point **&coordinateMesh,
                     double **&matrixPressure,
                     SystemPatemeters &systemParameters)
{
    const int MATRIX_PRESSURE_SIZE = systemParameters.n * systemParameters.n;
    const int MATRIX_CONTRIBUTION_SIZE = (systemParameters.n - 1) * (systemParameters.n - 1) * 2;

    initMatrix(matrixPressure, MATRIX_PRESSURE_SIZE, MATRIX_PRESSURE_SIZE);
    initMesh(coordinateMesh, systemParameters);

    initContributionMatrix(contributionMatrix, MATRIX_CONTRIBUTION_SIZE);

    createLocalMatrixForEveryElementHConst(contributionMatrix, coordinateMesh, systemParameters.n);
    createGlobalPressureMatrixHConst(matrixPressure, contributionMatrix, systemParameters.n);

    addBorderConditionsHConst(matrixPressure, systemParameters.n, MATRIX_PRESSURE_SIZE,
                              systemParameters.LOW_BORDER, systemParameters.HIGH_BORDER);

    outputPressureMatrix(matrixPressure, MATRIX_PRESSURE_SIZE);
}

/**
 * Create local stiffness matrixes for every element in consideration of its location. It can be top
 * and bottom triangle, so the formulas are different in these two cases.
 * */
void createLocalMatrixForEveryElementHConst(ContributionMatrix *&contributionMatrix,
                                            Point **&coordinateMesh,
                                            int n)
{
    int finiteElementNumber = -1;
    for (int i = 0; i < n - 1; i++)
    {
        for (int j = 0; j < n - 1; j++)
        {
            finiteElementNumber++;
            createLocalContributionMatrixForHConst(contributionMatrix[finiteElementNumber],
                                                   coordinateMesh[i][j],
                                                   coordinateMesh[i + 1][j],
                                                   coordinateMesh[i][j + 1]);
            finiteElementNumber++;
            createLocalContributionMatrixForHConst(contributionMatrix[finiteElementNumber],
                                                   coordinateMesh[i + 1][j + 1],
                                                   coordinateMesh[i][j + 1],
                                                   coordinateMesh[i + 1][j]);
        }
    }
}

/**
 * Fill local contribution matrix 
 * */
void createLocalContributionMatrixForHConst(ContributionMatrix localMatrix,
                                            Point pointI, Point pointJ, Point pointK)
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
            localMatrix.set(i, j, valB + valG);
        }
    }
}

/**
 * Count area of finite element
 * */
double countArea(Point pointI, Point pointJ, Point pointK)
{
    return fabs(0.5 * (pointJ.getX() * pointK.getY() - pointK.getX() * pointJ.getY() +
                       pointJ.getY() * pointI.getX() - pointK.getY() * pointI.getX() +
                       pointK.getX() * pointI.getY() - pointJ.getX() * pointI.getY()));
}

/**
 * Transform local matrixes and local indexes to global matrix and global indexes
 * */
void createGlobalPressureMatrixHConst(double **&matrixPressure, ContributionMatrix *&contributionMatrix, int n)
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

            for (int iterator1 = 0; iterator1 < 3; iterator1++)
                for (int iterator2 = 0; iterator2 < 3; iterator2++)
                    matrixPressure[globalNodeNumbersIJK[iterator1]][globalNodeNumbersIJK[iterator2]] +=
                        contributionMatrix[finiteElementNumber].get(iterator1, iterator2);

            finiteElementNumber++;

            //for bottom triangle
            globalNodeNumbersIJK[0] = globalNodeNumbersIJK[1] + 1;
            swap(globalNodeNumbersIJK[1], globalNodeNumbersIJK[2]);

            for (int iterator1 = 0; iterator1 < 3; iterator1++)
                for (int iterator2 = 0; iterator2 < 3; iterator2++)
                    matrixPressure[globalNodeNumbersIJK[iterator1]][globalNodeNumbersIJK[iterator2]] +=
                        contributionMatrix[finiteElementNumber].get(iterator1, iterator2);

            finiteElementNumber++;
        }
}

/**
 * The global stiffness matrix was created. In this function we take this matrix and change some
 * equations in it. For exaple, we have border condition on the right side of the region, in node n
 * (this is first node in secong row), we take equation number n, and just set value of node n to 
 * the condition(change the equation to x_n = p_n)
 * */
void addBorderConditionsHConst(double **&matrixResult,
                               int n,
                               int MATRIX_PRESSURE_SIZE,
                               int LOW_BORDER,
                               int HIGH_BORDER)
{
    double *rightPart = new double[MATRIX_PRESSURE_SIZE];

    //left
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < MATRIX_PRESSURE_SIZE; j++)
            matrixResult[i * n][j] = 0.0;

        matrixResult[i * n][i * n] = 1.0;
        rightPart[i * n] = LOW_BORDER;
    }

    //right
    for (int i = 1; i <= n; i++)
    {
        for (int j = 0; j < MATRIX_PRESSURE_SIZE; j++)
            matrixResult[i * n - 1][j] = 0.0;

        matrixResult[i * n - 1][i * n - 1] = 1.0;
        rightPart[i * n - 1] = LOW_BORDER;
    }

    //n row
    for (int i = MATRIX_PRESSURE_SIZE - n; i < MATRIX_PRESSURE_SIZE; i++)
    {
        for (int j = 0; j < MATRIX_PRESSURE_SIZE; j++)
            matrixResult[i][j] = 0.0;

        matrixResult[i][i] = 1.0;
        rightPart[i] = LOW_BORDER;
    }

    //0 row
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < MATRIX_PRESSURE_SIZE; j++)
            matrixResult[i][j] = 0.0;

        matrixResult[i][i] = 1.0;
        rightPart[i] = HIGH_BORDER;
    }

    fstream myFile;
    myFile.open(FILE_VECTOR_RIGHT_PART_OUTPUT, fstream::out);
    for (int i = 0; i < MATRIX_PRESSURE_SIZE; i++)
        myFile << rightPart[i] << endl;
}