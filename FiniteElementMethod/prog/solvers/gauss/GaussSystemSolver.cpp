#include "GaussSystemSolver.hpp"

const double epsilon = 1e-12;
std::string FILENAME_MATRIX_PRESSURE = "data/fem_output/pressureMatrix.txt";
std::string FILENAME_RIGHT_PART = "data/fem_output/rightPart.txt";
std::string FILENAME_SOLUTION = "data/gauss_output/solution.txt";

/**
 * Solve the given equation using Gauss method
 * */
void solveEquation(const int size)
{
    //Matrix of the system
    double **A;

    //Right part of the system
    double *B;

    //Vector of solution
    double *X;

    allocateMemory(A, B, X, size);

    int readSuccesfully = readData(FILENAME_MATRIX_PRESSURE, FILENAME_RIGHT_PART, A, B, size);

    if (readSuccesfully)
        if (solveWithGauss(A, B, X, size))
            WriteVector(FILENAME_SOLUTION, X, size);
        else
            std::cout << "Matrix A is degenerate " << std::endl;

    freeMemory(A, B, X, size);
}

/**
 * Allocate memory for all components of the method
 * */
void allocateMemory(double **&A,
                    double *&B,
                    double *&X,
                    const int &n)
{
    A = new double *[n];

    for (int i = 0; i < n; i++)
        A[i] = new double[n];

    B = new double[n];

    X = new double[n];
}

/**
 * Read matrix and right part from file, if everything is succesful return true
 * else return false
 * */
bool readData(const std::string fileNameMatrix,
              const std::string fileNameVector,
              double **&matrixA,
              double *&vectorB,
              const int &n)
{

    //Read matrix
    std::ifstream matrixFile;
    matrixFile.open(fileNameMatrix);

    if (!matrixFile.is_open())
    {
        std::cerr << "Error: file with matrix is not open" << std::endl;
        return false;
    }

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            matrixFile >> matrixA[i][j];

    matrixFile.close();

    //Read right part
    std::ifstream vectorFile;
    vectorFile.open(fileNameVector);

    if (!vectorFile.is_open())
    {
        std::cerr << "Error: file with vector is not open" << std::endl;
        return false;
    }

    for (int i = 0; i < n; i++)
        vectorFile >> vectorB[i];

    return true;
}

bool solveWithGauss(double **&A, double *&B, double *&X, const int &size)
{
    double m;
    double k;
    bool noProblems = true;

    //Only for testing puposes
    std::cout << "Target: " << size << std::endl;

    for (int i = 0; i < size; i++)
        X[i] = B[i];

    //Straight Gauss method
    for (int i = 0; i < size - 1; i++)
    {
        //Only for testing puposes
        if (i % 100 == 0)
            std::cout << "Current progress: " << i << std::endl;

        if (matrixIsPrepared(A, B, i, size))
        {
            m = A[i][i];
            B[i] /= m;

            //divide all row on diagonal element
            for (int j = i; j < size; j++)
                A[i][j] /= m;

            //beginning from element after diagonal substract elements of higher row from lower rows
            for (int i1 = i + 1; i1 < size; i1++)
            {
                k = A[i1][i];
                B[i1] -= k * B[i];

                for (int j1 = i; j1 < size; j1++)
                    A[i1][j1] -= k * A[i][j1];
            }
        }
        else
        {
            noProblems = false;
            break;
        }
    }

    if (matrixIsPrepared(A, B, size - 1, size) == false)
        noProblems = false;

    //Backwards Gauss method
    if (noProblems == true)
    {
        B[size - 1] /= A[size - 1][size - 1];
        A[size - 1][size - 1] /= A[size - 1][size - 1];

        //begin cycle till the first row
        for (int i = size - 1; i >= 0; i--)
            for (int j = i - 1; j >= 0; j--)
            {
                B[j] -= A[j][i] * B[i];
                A[j][i] = 0.0;
            }

        diagonalizeEquation(A, B, X, size);
        for (int i = 0; i < size; i++)
            std::swap(B[i], X[i]);
    }
    return noProblems;
}

/**
 * Check if current minor is ready for computation(does not have zero in position (0,0) of current minor)
 * */
bool matrixIsPrepared(double **&A,
                      double *&B,
                      const int &i,
                      const int &size)
{
    int indexOfMaximumInColumn = searchMaxInColumn(A, i, size);

    if (indexOfMaximumInColumn != i)
        swapRows(A, B, i, indexOfMaximumInColumn);

    if (isDegenerate(A, i, size) == false)
        return true;
    else
        return false;
}

/**
 * Place back rows and columns that we have swaped before
 * */
void diagonalizeEquation(double **&A,
                         double *&B,
                         double *&X,
                         const int &size)
{
    for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++)
            if (A[i][j] > epsilon && i != j)
                swapRows(A, B, i, j);
}

/**
 * Searching maximum element in the the first column of current minor
 * */
int searchMaxInColumn(double **&A,
                      const int &currentMinor,
                      const int &size)
{
    int maxI = currentMinor;
    double max = fabs(A[maxI][currentMinor]);
    for (int i = currentMinor; i < size; i++)
        if (fabs(A[i][currentMinor]) > max)
        {
            maxI = i;
            max = fabs(A[maxI][currentMinor]);
        }
    return maxI;
}

void swapRows(double **&A, double *&B, const int &i1, const int &i2)
{
    std::swap(A[i1], A[i2]);
    std::swap(B[i1], B[i2]);
}

/**
 * Check if matrix is degenerate
 * */
bool isDegenerate(double **&A,
                  const int &i,
                  const int &size)
{
    int numberOfZerosInRow = 0;
    for (int j = 0; j < size; j++)
        if (fabs(A[i][j]) < epsilon)
            numberOfZerosInRow++;

    if (numberOfZerosInRow == size)
        return true;

    return false;
}

/**
 * Delete allocated memory as we dont need it any more
 * */
void freeMemory(double **&A,
                double *&B,
                double *&X,
                const int &n)
{
    delete[] B;
    delete[] X;

    for (int i = 0; i < n; i++)
        delete[] A[i];
    delete[] A;
}

int WriteVector(std::string fileNameOutput,
                double *&vector,
                const int &n)
{
    std::ofstream fileOutput;
    fileOutput.open(fileNameOutput);

    for (int i = 0; i < n; i++)
        fileOutput << vector[i] << "\n";
    fileOutput << "\n";

    fileOutput.close();
    return 0;
}