#pragma once

/**
 * Data class for storing information about local contribution matrix
 * */
class SecondOrderRectangleContributionMatrix
{
private:
    const int ROW = 8;
    const int COLUMN = 8;
    double **matrix;

public:
    double** getMatrix()
    {
        return matrix;
    }

    double getElement(int i, int j)
    {
        return matrix[i][j];
    }

    void setElement(int i, int j, double value)
    {
        matrix[i][j] = value;
    }

    SecondOrderRectangleContributionMatrix()
    {
        matrix = new double *[ROW];
        for (int i = 0; i < ROW; i++)
            matrix[i] = new double[COLUMN];

        for (int i = 0; i < ROW; i++)
            for (int j = 0; j < COLUMN; j++)
                matrix[i][j] = 0.0;
    }
};