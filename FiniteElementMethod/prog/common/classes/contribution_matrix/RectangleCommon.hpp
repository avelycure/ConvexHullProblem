class RectnangleContributionMatrix
{
public:
    const int ROW = 4;
    const int COLUMN = 4;
    double **matrix;

    double getElement(int i, int j)
    {
        return matrix[i][j];
    }
    void setElement(int i, int j, double value)
    {
        matrix[i][j] = value;
    }
    RectnangleContributionMatrix()
    {
        matrix = new double *[ROW];
        for (int i = 0; i < ROW; i++)
            matrix[i] = new double[COLUMN];

        for (int i = 0; i < ROW; i++)
            for (int j = 0; j < COLUMN; j++)
                matrix[i][j] = 0.0;
    }
};

class RectnangleRightPart
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
    RectnangleRightPart()
    {
        vector = new double[ROWS];
        
        for (int i = 0; i < ROWS; i++)
            vector[i] = 0.0;
    }
};