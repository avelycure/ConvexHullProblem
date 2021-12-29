#pragma once

/**
 * This class represents point of the mesh
 * */
class Point
{
private:
    double x;
    double y;

public:
    Point(double xP, double yP)
    {
        x = xP;
        y = yP;
    }

    Point()
    {
        x = 0.0;
        y = 0.0;
    }

    double getX()
    {
        return x;
    }
    double getY()
    {
        return y;
    }
    void setX(double xParam)
    {
        x = xParam;
    }
    void setY(double yParam)
    {
        y = yParam;
    }
};