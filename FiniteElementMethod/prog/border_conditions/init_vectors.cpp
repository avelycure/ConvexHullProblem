#include <fstream>
#include <string>
#include <iostream>
#include <math.h>

using namespace std;

void outputVector(string fileName, double *vector, int size)
{
    fstream myFile;
    myFile.open(fileName, fstream::out);

    for (int i = 0; i < size; i++)
        myFile << vector[i] << endl;

    myFile.close();
}

void setValue(double *vector, int size, double value)
{
    for (int i = 0; i < size; i++)
        vector[i] = value;
}

void setValue(double *vector, int size)
{
    for (int i = 0; i < size; i++)
        vector[i] = sin(i) + rand() % 10 * 0.137;
}

int main()
{
    int size = 900; //n * n
    double *borderValues = new double[size];

    setValue(borderValues, size, 2.0);
    outputVector("top.txt", borderValues, size);

    setValue(borderValues, size, 1.0);
    outputVector("bottom.txt", borderValues, size);

    setValue(borderValues, size, 1.0);
    outputVector("right.txt", borderValues, size);

    setValue(borderValues, size, 1.0);
    outputVector("left.txt", borderValues, size);
    return 0;
}
