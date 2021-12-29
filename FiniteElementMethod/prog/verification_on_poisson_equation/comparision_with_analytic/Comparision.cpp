#include "Comparision.hpp"

double PI = 3.14159265358979323846;

/**
 * Function to compare solutions of all methods except rectangle second order
 * */
double compareWithAnalyticNormSum(Point **&coordinateMesh,
                                  int n,
                                  SystemParameters systemParameters)
{
    double norm = 0.0;
    double *solution = new double[n * n];

    readSolution("data/gauss_output/solution.txt", solution, n * n);

    for (int i = 1; i < n - 1; i++)
        for (int j = 1; j < n - 1; j++)
            norm += pow((solution[i * n + j] - analyticOfPoissonEquation(coordinateMesh[i][j].getX(),
                                                                         coordinateMesh[i][j].getY(),
                                                                         systemParameters)),
                        2.0);

    outputPoissonSolution(coordinateMesh,
                          "data/gauss_output/a_solution.txt",
                          n,
                          systemParameters);

    delete[] solution;
    return sqrt(norm);
}

double compareWithAnalyticNormMax(Point **&coordinateMesh,
                                  int n,
                                  SystemParameters systemParameters)
{
    double *solution = new double[n * n];
    readSolution("data/gauss_output/solution.txt", solution, n * n);
    double norm = 0.0;

    for (int i = 1; i < n - 1; i++)
        for (int j = 1; j < n - 1; j++)
            if (fabs(solution[i * n + j] - analyticOfPoissonEquation(coordinateMesh[i][j].getX(),
                                                                     coordinateMesh[i][j].getY(),
                                                                     systemParameters)) > norm)
            {
                norm = fabs(solution[i * n + j] - analyticOfPoissonEquation(coordinateMesh[i][j].getX(),
                                                                            coordinateMesh[i][j].getY(),
                                                                            systemParameters));
                //std::cout << i << " " << j << std::endl;
            }
    outputPoissonSolution(coordinateMesh,
                          "data/gauss_output/a_solution.txt",
                          n,
                          systemParameters);

    outputDifference(solution,
                     coordinateMesh,
                     "data/gauss_output/d_solution.txt",
                     n,
                     systemParameters);

    delete[] solution;
    return norm;
}

double compareWithAnalyticRectangleSecondOrderNormMax(Point **&coordinateMesh,
                                                      int n,
                                                      SystemParameters systemParameters)
{
    double norm = 0.0;
    int solutionLenght = systemParameters.n * systemParameters.n - ((systemParameters.n - 1) / 2) * ((systemParameters.n - 1) / 2);
    double *solution = new double[solutionLenght];

    readSolution("data/gauss_output/solution.txt", solution, solutionLenght);

    int node = n;

    for (int i = 1; i < n - 1; i++)
    {
        if (i % 2 == 0)
        {
            node++;
            for (int j = 1; j < n - 1; j++)
            {
                if (fabs(solution[node] - analyticOfPoissonEquation(coordinateMesh[i][j].getX(),
                                                                    coordinateMesh[i][j].getY(),
                                                                    systemParameters)) > norm)
                {
                    norm = fabs(solution[node] - analyticOfPoissonEquation(coordinateMesh[i][j].getX(),
                                                                           coordinateMesh[i][j].getY(),
                                                                           systemParameters));
                }
                //std::cout << "0Node: " << node << " (" << coordinateMesh[i][j].getX() << " " << coordinateMesh[i][j].getY() << ")" << std::endl;
                node++;
            }
        }

        if (i % 2 == 1)
        {
            node++;
            for (int j = 2; j < n - 2; j++)
            {
                if (j % 2 == 0)
                {

                    if (fabs(solution[node] - analyticOfPoissonEquation(coordinateMesh[i][j].getX(),
                                                                        coordinateMesh[i][j].getY(),
                                                                        systemParameters)) > norm)
                        norm = fabs(solution[node] - analyticOfPoissonEquation(coordinateMesh[i][j].getX(),
                                                                               coordinateMesh[i][j].getY(),
                                                                               systemParameters));

                    //std::cout << "1Node: " << node << " (" << coordinateMesh[i][j].getX() << " " << coordinateMesh[i][j].getY() << ")" << std::endl;
                    node++;
                }
            }
        }

        node++;
    }

    outputAnalyticSolutionOnSquareMesh(coordinateMesh,
                                       "data/gauss_output/a_solution.txt",
                                       n,
                                       systemParameters);

    outputDifferenceOnSquareMesh(solution,
                                 coordinateMesh,
                                 "data/gauss_output/d_solution.txt",
                                 n,
                                 systemParameters);

    return norm;
}

double compareWithAnalyticRectangleSecondOrderNormSum(Point **&coordinateMesh,
                                                      int n,
                                                      SystemParameters systemParameters)
{
    double norm = 0.0;
    int solutionLenght = systemParameters.n * systemParameters.n - ((systemParameters.n - 1) / 2) * ((systemParameters.n - 1) / 2);
    double *solution = new double[solutionLenght];

    readSolution("data/gauss_output/solution.txt", solution, solutionLenght);

    int node = n;

    for (int i = 1; i < n - 1; i++)
    {
        if (i % 2 == 0)
        {
            node++;
            for (int j = 1; j < n - 1; j++)
            {
                norm += pow(solution[node] - analyticOfPoissonEquation(coordinateMesh[i][j].getX(),
                                                                       coordinateMesh[i][j].getY(),
                                                                       systemParameters),
                            2.0);
                node++;
            }
        }

        if (i % 2 == 1)
        {
            node++;
            for (int j = 2; j < n - 2; j++)
            {
                if (j % 2 == 0)
                {
                    norm += pow(solution[node] - analyticOfPoissonEquation(coordinateMesh[i][j].getX(),
                                                                           coordinateMesh[i][j].getY(),
                                                                           systemParameters),
                                2.0);
                    node++;
                }
            }
        }

        node++;
    }

    outputAnalyticSolutionOnSquareMesh(coordinateMesh,
                                       "data/gauss_output/a_solution.txt",
                                       n,
                                       systemParameters);

    outputDifferenceOnSquareMesh(solution,
                                 coordinateMesh,
                                 "data/gauss_output/d_solution.txt",
                                 n,
                                 systemParameters);

    return sqrt(norm);
}

double analyticOfPoissonEquation(double x, double y, SystemParameters systemParameters)
{
    return systemParameters.poissonC1 * x * cos(systemParameters.poissonC2 * y);
}

/**
 * Read solution from file, if everything is succesful return true
 * else return false
 * */
bool readSolution(const std::string fileNameVector,
                  double *&solution,
                  const int &n)
{
    //Read solution
    std::ifstream solutionFile;
    solutionFile.open(fileNameVector);

    if (!solutionFile.is_open())
    {
        std::cerr << "Error: file with solution is not open" << std::endl;
        return false;
    }

    for (int i = 0; i < n; i++)
        solutionFile >> solution[i];

    return true;
}

void outputDifference(double *&solution,
                      Point **&coordinateMesh,
                      std::string fileNameOutput,
                      const int &n,
                      SystemParameters systemParameters)
{
    std::ofstream fileOutput;
    fileOutput.open(fileNameOutput);

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            fileOutput << fabs(solution[i * n + j] - analyticOfPoissonEquation(coordinateMesh[i][j].getX(),
                                                                               coordinateMesh[i][j].getY(),
                                                                               systemParameters))
                       << "\n";

    fileOutput.close();
}

void outputPoissonSolution(Point **&coordinateMesh,
                           std::string fileNameOutput,
                           const int &n,
                           SystemParameters systemParameters)
{
    std::ofstream fileOutput;
    fileOutput.open(fileNameOutput);

    int prec = std::numeric_limits<double>::digits10 + 2;
    int exponent_digits = std::log10(std::numeric_limits<double>::max_exponent10) + 1; // generally 3
    int exponent_sign = 1;                                                             // 1.e-123
    int exponent_symbol = 1;                                                           // 'e' 'E'
    int digits_sign = 1;
    int digits_dot = 1; // 1.2

    int division_extra_space = 1;
    int width = prec + exponent_digits + digits_sign + exponent_sign + digits_dot + exponent_symbol + division_extra_space;

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            fileOutput << std::setprecision(prec) << std::setw(width) << analyticOfPoissonEquation(coordinateMesh[i][j].getX(), coordinateMesh[i][j].getY(), systemParameters)
                       << "\n";

    fileOutput.close();
}

void outputAnalyticSolutionOnSquareMesh(Point **&coordinateMesh,
                                        std::string fileNameOutput,
                                        const int &n,
                                        SystemParameters systemParameters)
{
    std::ofstream fileOutput;
    fileOutput.open(fileNameOutput);

    int node = -1;

    for (int i = 0; i < n; i++)
    {
        if (i % 2 == 0)
        {
            node++;
            for (int j = 0; j < n; j++)
            {
                fileOutput << analyticOfPoissonEquation(coordinateMesh[i][j].getX(), coordinateMesh[i][j].getY(), systemParameters)
                           << "\n";
                node++;
            }
        }

        if (i % 2 == 1)
        {
            node++;
            for (int j = 0; j < n; j++)
            {
                if (j % 2 == 0)
                {
                    fileOutput << analyticOfPoissonEquation(coordinateMesh[i][j].getX(),
                                                            coordinateMesh[i][j].getY(),
                                                            systemParameters)
                               << "\n";
                    node++;
                }
            }
        }
        node++;
    }

    fileOutput.close();
}

void outputDifferenceOnSquareMesh(double *&solution,
                                  Point **&coordinateMesh,
                                  std::string fileNameOutput,
                                  const int &n,
                                  SystemParameters systemParameters)
{
    std::ofstream fileOutput;
    fileOutput.open(fileNameOutput);

    int node = 0;

    for (int i = 0; i < n; i++)
    {
        if (i % 2 == 0)
        {
            for (int j = 0; j < n; j++)
            {
                fileOutput << fabs(solution[node] - analyticOfPoissonEquation(coordinateMesh[i][j].getX(),
                                                                              coordinateMesh[i][j].getY(),
                                                                              systemParameters))
                           << "\n";
                node++;
            }
        }

        if (i % 2 == 1)
        {
            for (int j = 0; j < n; j++)
                if (j % 2 == 0)
                {
                    fileOutput << fabs(solution[node] - analyticOfPoissonEquation(coordinateMesh[i][j].getX(),
                                                                                  coordinateMesh[i][j].getY(),
                                                                                  systemParameters))
                               << "\n";
                    node++;
                }
        }
    }
    fileOutput.close();
}