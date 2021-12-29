#include "BorderConditions.hpp"

/**
 * Not for rectangle meshes
 * Check *
 * */
void addBorderConditionsOnPressureValues(SystemParameters systemParameters,
                                         double **&matrixResult,
                                         double *&rightPart,
                                         int n,
                                         int MATRIX_PRESSURE_SIZE)
{
    double h = systemParameters.borderLength / (n - 1);
    double x, y;

    //0 row
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < MATRIX_PRESSURE_SIZE; j++)
            matrixResult[i][j] = 0.0;

        matrixResult[i][i] = 1.0;

        x = 0.0;
        y = i * h;

        rightPart[i] = systemParameters.poissonC1 * x * cos(systemParameters.poissonC2 * y);
    }

    //left
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < MATRIX_PRESSURE_SIZE; j++)
            matrixResult[i * n][j] = 0.0;

        matrixResult[i * n][i * n] = 1.0;

        y = 0.0;
        x = i * h;

        rightPart[i * n] = systemParameters.poissonC1 * x * cos(systemParameters.poissonC2 * y);
    }

    //right
    for (int i = 1; i <= n; i++)
    {
        for (int j = 0; j < MATRIX_PRESSURE_SIZE; j++)
            matrixResult[i * n - 1][j] = 0.0;

        matrixResult[i * n - 1][i * n - 1] = 1.0;

        //1.0
        y = systemParameters.borderLength;
        x = (i - 1) * h;

        rightPart[i * n - 1] = systemParameters.poissonC1 * x * cos(systemParameters.poissonC2 * y);
    }

    //n row
    int yi = 0;
    for (int i = MATRIX_PRESSURE_SIZE - n; i < MATRIX_PRESSURE_SIZE; i++)
    {
        for (int j = 0; j < MATRIX_PRESSURE_SIZE; j++)
            matrixResult[i][j] = 0.0;

        matrixResult[i][i] = 1.0;

        //1.0
        x = systemParameters.borderLength;
        y = yi * h;
        yi++;

        rightPart[i] = systemParameters.poissonC1 * x * cos(systemParameters.poissonC2 * y);
    }

    int prec = std::numeric_limits<double>::digits10 + 2;
    int exponent_digits = std::log10(std::numeric_limits<double>::max_exponent10) + 1; // generally 3
    int exponent_sign = 1;                                                             // 1.e-123
    int exponent_symbol = 1;                                                           // 'e' 'E'
    int digits_sign = 1;
    int digits_dot = 1; // 1.2

    int division_extra_space = 1;
    int width = prec + exponent_digits + digits_sign + exponent_sign + digits_dot + exponent_symbol + division_extra_space;

    std::fstream myFile1;
    myFile1.open("data/fem_output/rightPart.txt", std::fstream::out);
    for (int i = 0; i < MATRIX_PRESSURE_SIZE; i++)
        myFile1 << std::setprecision(prec) << std::setw(width) << rightPart[i] << std::endl;

    myFile1.close();
}

/**
 * For rectangle meshes
 * Fix yi like in previous + another
 * */
void addBorderConditionsForRectangleSecondOrder(double **&matrixResult,
                                                double *&rightPart,
                                                int n,
                                                int MATRIX_PRESSURE_SIZE,
                                                SystemParameters systemParameters)
{
    int nSkippedRows = 0;
    int numberOfFE = ((n - 1) / 2);

    double h = systemParameters.borderLength / (n - 1);
    double x, y;

    //left
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < MATRIX_PRESSURE_SIZE; j++)
            matrixResult[i * n - nSkippedRows * numberOfFE][j] = 0.0;

        matrixResult[i * n - nSkippedRows * numberOfFE][i * n - nSkippedRows * numberOfFE] = 1.0;

        y = 0.0;
        x = i * h;

        rightPart[i * n - nSkippedRows * numberOfFE] = systemParameters.poissonC1 * x * cos(systemParameters.poissonC2 * y);

        if (i % 2 == 1)
            nSkippedRows++;
    }

    nSkippedRows = 0;
    //right
    for (int i = 1; i <= n; i++)
    {
        //increase begining from second row
        if (i % 2 == 0)
            nSkippedRows++;

        for (int j = 0; j < MATRIX_PRESSURE_SIZE; j++)
            matrixResult[i * n - nSkippedRows * numberOfFE - 1][j] = 0.0;

        matrixResult[i * n - nSkippedRows * numberOfFE - 1][i * n - nSkippedRows * numberOfFE - 1] = 1.0;

        y = 1.0;
        x = (i - 1) * h;

        rightPart[i * n - nSkippedRows * numberOfFE - 1] = systemParameters.poissonC1 * x * cos(systemParameters.poissonC2 * y);
    }

    //n row
    int yi = 0;
    for (int i = MATRIX_PRESSURE_SIZE - n; i < MATRIX_PRESSURE_SIZE; i++)
    {
        for (int j = 0; j < MATRIX_PRESSURE_SIZE; j++)
            matrixResult[i][j] = 0.0;

        y = yi * h;
        x = 1.0;
        yi++;

        matrixResult[i][i] = 1.0;
        rightPart[i] = systemParameters.poissonC1 * x * cos(systemParameters.poissonC2 * y);
    }

    //0 row
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < MATRIX_PRESSURE_SIZE; j++)
            matrixResult[i][j] = 0.0;

        y = i * h;
        x = 0.0;

        matrixResult[i][i] = 1.0;
        rightPart[i] = systemParameters.poissonC1 * x * cos(systemParameters.poissonC2 * y);
    }

    std::fstream myFile1;
    myFile1.open("data/fem_output/rightPart.txt", std::fstream::out);
    for (int i = 0; i < MATRIX_PRESSURE_SIZE; i++)
        myFile1 << rightPart[i] << std::endl;

    myFile1.close();
}