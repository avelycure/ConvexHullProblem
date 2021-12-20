#pragma once

/**
 * Right part of the local contribution matrix using rectangle finite elements
 * */
class SecondOrderRectangleRightPart
{
private:
    const int ROWS = 8;
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

    SecondOrderRectangleRightPart()
    {
        vector = new double[ROWS];

        for (int i = 0; i < ROWS; i++)
            vector[i] = 0.0;
    }
};