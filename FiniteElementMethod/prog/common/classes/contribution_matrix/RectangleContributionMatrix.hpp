class RectangleContributionMatrix
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
    RectangleContributionMatrix()
    {
        matrix = new double *[ROW];
        for (int i = 0; i < ROW; i++)
            matrix[i] = new double[COLUMN];

        for (int i = 0; i < ROW; i++)
            for (int j = 0; j < COLUMN; j++)
                matrix[i][j] = 0.0;
    }
};