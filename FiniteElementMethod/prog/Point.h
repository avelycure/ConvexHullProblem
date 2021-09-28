class Point
{
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
    void setX(double xP)
    {
        x = xP;
    }
    void setY(double yP)
    {
        y = yP;
    }
};