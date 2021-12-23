#include "Comparision.hpp"

double compareWithAnalyticNormSum(double *&solutionFE, int n, SystemParameters systemParameters)
{
    double norm = 0;

    for (int i = 0; i < n; i++)
        norm += fabs(solutionFE[i] - analytic(1, 1));

    return norm;
}

double compareWithAnalyticNormMax(double *&solutionFE, int n, SystemParameters systemParameters)
{
    double norm = 0.0;

    for (int i = 0; i < n; i++)
        if (fabs(solutionFE[i] - analytic(1, 1)) > norm)
            norm = fabs(solutionFE[i] - analytic(1, 1));

    return norm;
}

double analytic(double x, double y)
{
    return x + y;
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