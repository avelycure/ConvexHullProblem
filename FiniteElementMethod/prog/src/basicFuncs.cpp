#include "../header.hpp"

void initMatrix(double **&matrix, int row, int column)
{
    matrix = new double *[row];
    for (int i = 0; i < row; i++)
        matrix[i] = new double[column];
}

void displayMatrix(double **matrix, int row, int column)
{
    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < column; j++)
            cout << matrix[i][j] << " ";

        cout << endl;
    }
}

void displayMesh(Point **coordinateMesh, int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
            cout << "(" << coordinateMesh[i][j].getX() << "," << coordinateMesh[i][j].getY() << ") ";

        cout << endl;
    }
}

void displayVector(double *mVector, int n)
{
    for (int i = 0; i < n; i++)
        cout << mVector[i] << " ";
    cout << endl;
}

/**
 * Display matrixes of all finit elements
 * */
void displayAllLocalMatrixes(ContributionMatrix *&ContributionMatrixParam, int n)
{
    for (int i = 0; i < n; i++)
    {
        cout << "matrix[" << i << "]" << endl;
        for (int ii = 0; ii < 3; ii++)
        {
            for (int jj = 0; jj < 3; jj++)
                cout << ContributionMatrixParam[i].get(ii, jj) << " ";

            cout << endl;
        }
        cout << endl;
    }
}

void outputPressureMatrix(double **matrixPressure, int MATRIX_PRESSURE_SIZE)
{
    fstream myFile;

    myFile.open(FILE_OUTPUT_MATRIX, fstream::out);
    for (int i = 0; i < MATRIX_PRESSURE_SIZE; i++)
    {
        for (int j = 0; j < MATRIX_PRESSURE_SIZE; j++)
            myFile << matrixPressure[i][j] << " ";

        myFile << endl;
    }
}

void outputVector(string fileName, double *vector, int size)
{
    fstream myFile;
    myFile.open(fileName, fstream::out);

    for (int i = 0; i < size; i++)
        myFile << vector[i] << endl;

    myFile.close();
}

void initVector(double *&p, int n)
{
    p = new double[n];
}

void initContributionMatrix(ContributionMatrix *&contributionMatrix, int MATRIX_CONTRIBUTION_SIZE)
{
    contributionMatrix = new ContributionMatrix[MATRIX_CONTRIBUTION_SIZE];
}

void initRightPart(RightPart *&localRigthParts, int MATRIX_CONTRIBUTION_SIZE)
{
    localRigthParts = new RightPart[MATRIX_CONTRIBUTION_SIZE];
}

/**
 * Set coordinates of the nodes
 * */
void initMesh(Point **&coordinateMesh, SystemPatemeters &systemParameters)
{
    int n = systemParameters.n;
    double xOrigin = systemParameters.xOrigin;
    double yOrigin = systemParameters.yOrigin;
    double h = systemParameters.L / (n - 1);
    coordinateMesh = new Point *[n];
    for (int i = 0; i < n; i++)
        coordinateMesh[i] = new Point[n];

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
        {
            coordinateMesh[i][j].setX(xOrigin + i * h);
            coordinateMesh[i][j].setY(yOrigin + j * h);
        }
}

/**
 * Read parameters of the system from json file. This done not to recompile program every time
 * we need to change parameters
 * */
void readSystemParameters(SystemPatemeters &systemParameters, string &method)
{
    nlohmann::json j;
    fstream fileInputSystem;
    fileInputSystem.open(FILE_SYSTEM_NAME);
    fileInputSystem >> j;
    int k = j["system"];
    fileInputSystem.close();

    nlohmann::json jsonObject;
    fstream fileInput;
    fileInput.open(FILE_PARAMETERS_NAME);
    fileInput >> jsonObject;

    method = jsonObject["configs"][k]["method"];
    systemParameters.n = jsonObject["configs"][k]["n"];
    systemParameters.L = jsonObject["configs"][k]["L"];
    systemParameters.U = jsonObject["configs"][k]["U"];
    systemParameters.mu = jsonObject["configs"][k]["mu"];
    systemParameters.xOrigin = jsonObject["configs"][k]["xOrigin"];
    systemParameters.yOrigin = jsonObject["configs"][k]["yOrigin"];
    systemParameters.borderConditions = jsonObject["configs"][k]["borderConditions"];

    systemParameters.hMin = jsonObject["configs"][k]["realInputGap"];
    systemParameters.LOW_BORDER = jsonObject["configs"][k]["lowBorder"];
    systemParameters.HIGH_BORDER = jsonObject["configs"][k]["highBorder"];
    systemParameters.k = jsonObject["configs"][k]["realK"];

    systemParameters.realPressure = jsonObject["configs"][k]["realPressure"];
    systemParameters.realK = jsonObject["configs"][k]["realK"];
    systemParameters.realInputGap = jsonObject["configs"][k]["realInputGap"];
    fileInput.close();
}

/**
 * This function transforms input parameters into dimensionless, so computation will be more
 * efficient, and we could see different models
 * */
void dimensionlessSystemParameters(SystemPatemeters &systemParameters, string &method)
{
    systemParameters.hMin /= systemParameters.realInputGap;
    systemParameters.LOW_BORDER /= systemParameters.realPressure;
    systemParameters.HIGH_BORDER /= systemParameters.realPressure;

    if (method == H_LINEAR)
        systemParameters.k = (systemParameters.k * systemParameters.L) / systemParameters.realInputGap;
    else
        systemParameters.k = 0.0;
}

/**
 * The global stiffness matrix was created. In this function we take this matrix and change some
 * equations in it. For exaple, we have border condition on the right side of the region, in node n
 * (this is first node in secong row), we take equation number n, and just set value of node n to 
 * the condition(change the equation to x_n = p_n)
 * */
void addBorderConditions(double **&matrixResult,
                         int n,
                         int MATRIX_PRESSURE_SIZE,
                         double LOW_BORDER,
                         double HIGH_BORDER)
{
    double *rightPart = new double[MATRIX_PRESSURE_SIZE];
    for (int i = 0; i < MATRIX_PRESSURE_SIZE; i++)
        rightPart[i] = 0.0;

    //0 row
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < MATRIX_PRESSURE_SIZE; j++)
            matrixResult[i][j] = 0.0;

        matrixResult[i][i] = 1.0;
        rightPart[i] = HIGH_BORDER;
    }

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

    fstream myFile;
    myFile.open(FILE_VECTOR_RIGHT_PART_OUTPUT, fstream::out);
    for (int i = 0; i < MATRIX_PRESSURE_SIZE; i++)
        myFile << rightPart[i] << endl;
}

/**
 * Overloaded function. Supports users input vector of border conditions
 * */
void addBorderConditions(double **&matrixResult,
                         double *&borderValues,
                         int n,
                         int MATRIX_PRESSURE_SIZE,
                         string borderPosition)
{
    /**
     * Vector of right parts
     * */
    double *rightPart = new double[MATRIX_PRESSURE_SIZE];

    if (borderPosition == "TOP")
    {
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < MATRIX_PRESSURE_SIZE; j++)
                matrixResult[i][j] = 0.0;

            matrixResult[i][i] = 1.0;
            rightPart[i] = borderValues[i];
        }
    }

    if (borderPosition == "BOTTOM")
    {
        for (int i = MATRIX_PRESSURE_SIZE - n; i < MATRIX_PRESSURE_SIZE; i++)
        {
            for (int j = 0; j < MATRIX_PRESSURE_SIZE; j++)
                matrixResult[i][j] = 0.0;

            matrixResult[i][i] = 1.0;
            rightPart[i] = borderValues[i];
        }
    }

    if (borderPosition == "RIGHT")
    {
        for (int i = 1; i <= n; i++)
        {
            for (int j = 0; j < MATRIX_PRESSURE_SIZE; j++)
                matrixResult[i * n - 1][j] = 0.0;

            matrixResult[i * n - 1][i * n - 1] = 1.0;
            rightPart[i * n - 1] = borderValues[i];
        }
    }

    if (borderPosition == "LEFT")
    {
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < MATRIX_PRESSURE_SIZE; j++)
                matrixResult[i * n][j] = 0.0;

            matrixResult[i * n][i * n] = 1.0;
            rightPart[i * n] = borderValues[i];
        }
    }

    fstream myFile;
    myFile.open(FILE_VECTOR_RIGHT_PART_OUTPUT, fstream::out);
    for (int i = 0; i < MATRIX_PRESSURE_SIZE; i++)
        myFile << rightPart[i] << endl;
}

void inputVector(const string fileNameVector, double *&borderValues, const int n)
{
    ifstream vectorFile;
    vectorFile.open(fileNameVector);

    if (!vectorFile.is_open())
        std::cerr << "Error: file with vector is not open\n";

    for (int i = 0; i < n; ++i)
        vectorFile >> borderValues[i];

    vectorFile.close();
}