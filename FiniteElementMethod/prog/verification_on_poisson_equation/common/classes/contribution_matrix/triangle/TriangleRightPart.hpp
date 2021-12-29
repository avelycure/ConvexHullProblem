#pragma once

/**
 * Right part of the local contribution matrix using triangle finite elements
 * */
class TriangleRightPart
{
private:
    const int COLUMN = 3;
    double *vector;

public:
    double getElement(int i)
    {
        return vector[i];
    }

    void setElement(int i, double value)
    {
        vector[i] = value;
    }

    TriangleRightPart()
    {
        vector = new double[COLUMN];
    }
};