class ContributionMatrix
{
public:
    const int ROW = 3;
    const int COLUMN = 3;
    double **matrix;

    double getElement(int i, int j)
    {
        return matrix[i][j];
    }
    void setElement(int i, int j, double value)
    {
        matrix[i][j] = value;
    }
    ContributionMatrix()
    {
        matrix = new double *[ROW];
        for (int i = 0; i < ROW; i++)
            matrix[i] = new double[COLUMN];
    }
};

class RightPart
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
    RightPart()
    {
        vector = new double[COLUMN];
    }
};