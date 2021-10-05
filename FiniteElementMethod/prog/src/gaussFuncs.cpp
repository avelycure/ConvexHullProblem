#include "../header.hpp"

const double epsilon = 1e-12;
string FILENAME_MATRIX_PRESSURE = "data/pressureMatrix.txt";
string FILENAME_RIGHT_PART = "data/rightPart.txt";
string FILENAME_SOLUTION = "result/solution.txt";

int solveEquation(const int size)
{
    double **A;
    double *B;
    double *X;

    allocateMemory(A, B, X, size);
    readData(FILENAME_MATRIX_PRESSURE, FILENAME_RIGHT_PART, A, B, size);

    if (gaussMethod(A, B, X, size))
        writeVector(FILENAME_SOLUTION, X, size);
    else
        cout << "Matrix A is degenerate " << endl;

    freeMemory(A, B, X, size);
    return 0;
}

int allocateMemory(double **&A, double *&B, double *&X, const int &n)
{
    A = new double *[n];

    for (int i = 0; i < n; ++i)
        A[i] = new double[n];

    B = new double[n];
    X = new double[n];
    return 0;
}

int readData(const string fileNameMatrix, const string fileNameVector, double **&matrixA, double *&vectorB, const int &n)
{
    ifstream matrixFile;
    matrixFile.open(fileNameMatrix);

    if (!matrixFile.is_open())
    {
        cerr << "Error: file with matrix is not open\n";
        return 1;
    }

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
            matrixFile >> matrixA[i][j];
    }

    matrixFile.close();

    ifstream vectorFile;
    vectorFile.open(fileNameVector);

    if (!vectorFile.is_open())
    {
        cerr << "Error: file with vector is not open\n";
        return 1;
    }

    for (int i = 0; i < n; ++i)
        vectorFile >> vectorB[i];

    return 0;
}

bool gaussMethod(double **&A, double *&B, double *&X, const int &size)
{
    double m;
    double k;
    bool noProblems = true;
    vector<tuple<int, int>> permutations = {};

    cout << "Target: " << size << endl;

    for (int i = 0; i < size; ++i)
        X[i] = B[i];

    for (int i = 0; i < size - 1; i++)
    {

        if (i % 100 == 0)
            cout << "Current progress: "<< i << endl;

        if (matrixIsPrepared(A, B, i, permutations, size))
        {
            m = A[i][i];
            B[i] /= m;

            for (int j = i; j < size; j++)
                A[i][j] /= m;

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

    if (matrixIsPrepared(A, B, size - 1, permutations, size) == false)
        noProblems = false;

    if (noProblems == true)
    {
        B[size - 1] /= A[size - 1][size - 1];
        A[size - 1][size - 1] /= A[size - 1][size - 1];
        for (int i = size - 1; i > -1; i--)
            for (int j = i - 1; j > -1; j--)
            {
                B[j] -= A[j][i] * B[i];
                A[j][i] = 0.0;
            }

        diagonalizeEquation(A, B, X, size, permutations);
        for (int i = 0; i < size; i++)
            swap(B[i], X[i]);
    }
    permutations.clear();
    return noProblems;
}

bool matrixIsPrepared(double **&A, double *&B, const int &i, vector<tuple<int, int>> &permutations, const int &size)
{
    tuple<int, int> t;
    t = searchMax(A, i, size);

    if (get<0>(t) != i)
        swapRows(A, B, i, get<0>(t));
    if (isDegenerate(A, i, size) == false)
        return true;
    else
        return false;
}

void diagonalizeEquation(double **&A, double *&B, double *&X, const int &size, vector<tuple<int, int>> &permutations)
{
    for (int i = permutations.size() - 1; i > -1; i--)
        swapColomns(A, get<0>(permutations[i]), get<1>(permutations[i]), size);

    for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++)
            if (A[i][j] > epsilon && i != j)
                swapRows(A, B, i, j);
}

tuple<int, int> searchMax(double **&A, const int &currentMinor, const int &size)
{
    int max_i = currentMinor;
    double max = fabs(A[max_i][currentMinor]);
    for (int i = currentMinor; i < size; i++)
        if (fabs(A[i][currentMinor]) > max)
        {
            max_i = i;
            max = fabs(A[max_i][currentMinor]);
        }
    return make_tuple(max_i, currentMinor);
}

void swapRows(double **&A, double *&B, const int &i1, const int &i2)
{
    swap(A[i1], A[i2]);
    swap(B[i1], B[i2]);
}

void swapColomns(double **&A, const int &j1, const int &j2, const int &size)
{
    for (int i = 0; i < size; i++)
        swap(A[i][j1], A[i][j2]);
}

bool isDegenerate(double **&A, const int &i, const int &size)
{
    int k = 0;
    for (int j = 0; j < size; j++)
        if (fabs(A[i][j]) < epsilon)
            k++;

    if (k == size)
        return true;

    return false;
}

void freeMemory(double **&A, double *&B, double *&X, const int &n)
{
    delete[] B;
    delete[] X;

    for (int i = 0; i < n; ++i)
        delete[] A[i];
    delete[] A;
}

void writeVector(string fileNameOutput, double *&vector, const int &n)
{
    ofstream fileOutput;
    fileOutput.open(fileNameOutput);

    for (int i = 0; i < n; ++i)
        fileOutput << vector[i] << "\n";
    fileOutput << "\n";

    fileOutput.close();
}