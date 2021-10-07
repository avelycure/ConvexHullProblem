/**
 * Class represents stiffness matrix
 * */
class ContributionMatrix
{
    double **matrix;

public:
    const int ROW = 3;
    const int COLUMN = 3;

    double get(int i, int j)
    {
        return matrix[i][j];
    }

    void set(int i, int j, double value)
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

/**
 * Class represents local right part
 * */
class RightPart
{
    double *vector;

public:
    const int COLUMN = 3;

    double get(int i)
    {
        return vector[i];
    }
    void set(int i, double value)
    {
        vector[i] = value;
    }
    RightPart()
    {
        vector = new double[COLUMN];
    }
};