class TriangleRightPart
{
public:
    const int COLUMN = 3;
    double *vector;

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