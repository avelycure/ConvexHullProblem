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
            norm += pow(solution[i * n + j] - analyticOfLaplasEquation(coordinateMesh[i][j].getX(),
                                                                        coordinateMesh[i][j].getY(),
                                                                        systemParameters),2.0);

    outputLaplasSolution(coordinateMesh,
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
            if (fabs(solution[i * n + j] - analyticOfLaplasEquation(coordinateMesh[i][j].getX(),
                                                                    coordinateMesh[i][j].getY(),
                                                                    systemParameters)) > norm)
            {
                norm = fabs(solution[i * n + j] - analyticOfLaplasEquation(coordinateMesh[i][j].getX(),
                                                                           coordinateMesh[i][j].getY(),
                                                                           systemParameters));
                //std::cout << i << " " << j << std::endl;
            }
    outputLaplasSolution(coordinateMesh,
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
                if (fabs(solution[node] - analyticOfLaplasEquation(coordinateMesh[i][j].getX(),
                                                                   coordinateMesh[i][j].getY(),
                                                                   systemParameters)) > norm)
                {
                    norm = fabs(solution[node] - analyticOfLaplasEquation(coordinateMesh[i][j].getX(),
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

                    if (fabs(solution[node] - analyticOfLaplasEquation(coordinateMesh[i][j].getX(),
                                                                       coordinateMesh[i][j].getY(),
                                                                       systemParameters)) > norm)
                        norm = fabs(solution[node] - analyticOfLaplasEquation(coordinateMesh[i][j].getX(),
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
                norm += pow(solution[node] - analyticOfLaplasEquation(coordinateMesh[i][j].getX(),
                                                                       coordinateMesh[i][j].getY(),
                                                                       systemParameters),2.0);
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
                    norm += pow(solution[node] - analyticOfLaplasEquation(coordinateMesh[i][j].getX(),
                                                                           coordinateMesh[i][j].getY(),
                                                                           systemParameters),2.0);
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
    
    delete[] solution;
    return sqrt(norm);
}

double analyticOfLaplasEquation(double x, double y, SystemParameters systemParameters)
{
    //x = 0
    double f1 = systemParameters.HIGH_BORDER;
    //x = 1
    double f2 = systemParameters.LOW_BORDER;
    //y = 0
    double f3 = systemParameters.LOW_BORDER;
    //y = 1
    double f4 = systemParameters.LOW_BORDER;

    double sum1 = 0.0;
    double sum2 = 0.0;
    double sum3 = 0.0;
    double sum4 = 0.0;

    int N = 100;

    for (int n = 1; n <= N; n++)
        sum1 += a(n) * sinh(n * PI * (1.0 - x)) * sin(n * PI * y);

    for (int n = 1; n <= N; n++)
        sum2 += b(n) * sinh(n * PI * x) * sin(n * PI * y);

    for (int n = 1; n <= N; n++)
        sum3 += c(n) * sin(n * PI * x) * sinh(n * PI * (1.0 - y));

    for (int n = 1; n <= N; n++)
        sum4 += d(n) * sin(n * PI * x) * sinh(n * PI * y);

    return sum1 + sum2 + sum3 + sum4;
}

double a(int n)
{
    return (2.0 * (2.0 - 2.0 * cos(n * PI))) / (n * PI * sinh(n * PI));
}

double b(int n)
{
    return (2.0 * (1.0 - cos(n * PI))) / (n * PI * sinh(n * PI));
}

double c(int n)
{
    return (2.0 * (1.0 - cos(n * PI))) / (n * PI * sinh(n * PI));
}

double d(int n)
{
    return (2.0 * (1.0 - cos(n * PI))) / (n * PI * sinh(n * PI));
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
            if (i == 0 || j == 0 || i == n - 1 || j == n - 1)
                fileOutput << 0.0 << "\n";
            else
                fileOutput << fabs(solution[i * n + j] - analyticOfLaplasEquation(coordinateMesh[i][j].getX(),
                                                                                  coordinateMesh[i][j].getY(),
                                                                                  systemParameters))
                           << "\n";

    fileOutput.close();
}

void outputLaplasSolution(Point **&coordinateMesh,
                          std::string fileNameOutput,
                          const int &n,
                          SystemParameters systemParameters)
{
    std::ofstream fileOutput;
    fileOutput.open(fileNameOutput);

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            if (i == 0)
                fileOutput << systemParameters.HIGH_BORDER << "\n";
            else if (j == 0 && i != 0)
                fileOutput << systemParameters.LOW_BORDER << "\n";
            else if (j == n - 1 && i != 0)
                fileOutput << systemParameters.LOW_BORDER << "\n";
            else if (i == n - 1)
                fileOutput << systemParameters.LOW_BORDER << "\n";
            else
                fileOutput << analyticOfLaplasEquation(coordinateMesh[i][j].getX(),
                                                       coordinateMesh[i][j].getY(),
                                                       systemParameters)
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
                if (i == 0 && j == 0)
                    fileOutput << 0.0 << "\n";
                else if (i == n - 1 && j == 0)
                    fileOutput << 0.0 << "\n";
                else if (i == n - 1 && j == n - 1)
                    fileOutput << 0.0 << "\n";
                else if (i == 0 && j == n - 1)
                    fileOutput << 0.0 << "\n";
                else
                    fileOutput << analyticOfLaplasEquation(coordinateMesh[i][j].getX(), coordinateMesh[i][j].getY(), systemParameters)
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
                    if (i == 0 && j == 0)
                        fileOutput << 0.0 << "\n";
                    else if (i == n - 1 && j == 0)
                        fileOutput << 0.0 << "\n";
                    else if (i == n - 1 && j == n - 1)
                        fileOutput << 0.0 << "\n";
                    else if (i == 0 && j == n - 1)
                        fileOutput << 0.0 << "\n";
                    else
                        fileOutput << analyticOfLaplasEquation(coordinateMesh[i][j].getX(),
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
                if (i == 0 || j == 0 || i == n - 1 || j == n - 1)
                    fileOutput << 0.0 << "\n";
                else
                    fileOutput << fabs(solution[node] - analyticOfLaplasEquation(coordinateMesh[i][j].getX(),
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
                    if (i == 0 || j == 0 || i == n - 1 || j == n - 1)
                        fileOutput << 0.0 << "\n";
                    else
                        fileOutput << fabs(solution[node] - analyticOfLaplasEquation(coordinateMesh[i][j].getX(),
                                                                                     coordinateMesh[i][j].getY(),
                                                                                     systemParameters))
                                   << "\n";
                    node++;
                }
        }
    }
    fileOutput.close();
}