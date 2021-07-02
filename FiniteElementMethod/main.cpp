#include "header.hpp"
int main()
{
    string method;
    double *rightPart;
    Point **coordinateMesh;
    double **matrixPressure;
    RightPart *localRigthParts;
    SystemPatemeters systemParameters;
    ContributionMatrix *contributionMatrix;

    readSystemParameters(systemParameters, method);

    if (method == H_CONST)
        solveWithHConst(contributionMatrix, coordinateMesh, matrixPressure,
                        systemParameters);

    if (method == H_LINEAR)
        solveWithHLinear(contributionMatrix, localRigthParts, coordinateMesh,
                         matrixPressure, rightPart, systemParameters);
    solveEquation(systemParameters.n * systemParameters.n);
}
