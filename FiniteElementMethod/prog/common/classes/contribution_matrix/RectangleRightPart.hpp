class RectangleRightPart
{
public:
    const int ROWS = 4;
    double *vector;

    double getElement(int i)
    {
        return vector[i];
    }
    void setElement(int i, double value)
    {
        vector[i] = value;
    }
    RectangleRightPart()
    {
        vector = new double[ROWS];
        
        for (int i = 0; i < ROWS; i++)
            vector[i] = 0.0;
    }
};