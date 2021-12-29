#pragma once

/**
 * 
 * */
class TriangleContributionMatrixSecondOrder
{
public:
    const int ROW = 6;
    const int COLUMN = 6;
    double **matrix;

    double getElement(int i, int j)
    {
        return matrix[i][j];
    }
    void setElement(int i, int j, double value)
    {
        matrix[i][j] = value;
    }
    TriangleContributionMatrixSecondOrder()
    {
        matrix = new double *[ROW];
        for (int i = 0; i < ROW; i++)
            matrix[i] = new double[COLUMN];
    }
};