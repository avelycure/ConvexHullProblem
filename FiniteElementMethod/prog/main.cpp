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

    if (method == "der")
        solveWithHLinearWithDerBC(contributionMatrix, localRigthParts, coordinateMesh,
                                  matrixPressure, rightPart, systemParameters);

    if(method == "derh")
        solveWithHConstBCLR(contributionMatrix, coordinateMesh, matrixPressure,
                        systemParameters);

    solveEquation(systemParameters.n * systemParameters.n);
}
