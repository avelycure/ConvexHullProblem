#include "BorderConditions.hpp"

void addBorderConditionsOnPressureValues(double **&matrixResult,
                                         double *&rightPart,
                                         int n,
                                         int MATRIX_PRESSURE_SIZE,
                                         double LOW_BORDER,
                                         double HIGH_BORDER)
{
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

void addBorderConditionsForRectangleSecondOrder(double **&matrixResult,
                                                double *&rightPart,
                                                int n,
                                                int MATRIX_PRESSURE_SIZE,
                                                double LOW_BORDER,
                                                double HIGH_BORDER)
{
    int nSkippedRows = 0;
    int numberOfFE = ((n - 1) / 2);
    //left
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < MATRIX_PRESSURE_SIZE; j++)
            matrixResult[i * n - nSkippedRows * numberOfFE][j] = 0.0;

        matrixResult[i * n - nSkippedRows * numberOfFE][i * n - nSkippedRows * numberOfFE] = 1.0;
        rightPart[i * n - nSkippedRows * numberOfFE] = LOW_BORDER;

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
        rightPart[i * n - nSkippedRows * numberOfFE - 1] = LOW_BORDER;
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

    std::fstream myFile1;
    myFile1.open("data/fem_output/rightPart.txt", std::fstream::out);
    for (int i = 0; i < MATRIX_PRESSURE_SIZE; i++)
        myFile1 << rightPart[i] << std::endl;

    myFile1.close();
}