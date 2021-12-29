#pragma once
/**
 * Right part of the local contribution matrix using triangle finite elements
 * */
class TriangleRightPartSecondOrder
{
private:
    const int size = 6;
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

    TriangleRightPartSecondOrder()
    {
        vector = new double[size];
    }
};