#include "main.hpp"
int main()
{
    string method;
    double *rightPart;
    Point **coordinateMesh;
    double **matrixPressure;
    RightPart *localRigthParts;

    SystemPatemeters systemParameters;
    ContributionMatrix *contributionMatrix;

    RectnangleContributionMatrix *contributionMatrixR;
    RectnangleRightPart *localRigthPartsR;

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

    if (method == "derh")
        solveWithHConstBCLR(contributionMatrix, coordinateMesh, matrixPressure,
                            systemParameters);

    if (method == "rect")
        solveWithRectangleFiniteElements(contributionMatrixR, localRigthPartsR, coordinateMesh,
                                         matrixPressure, rightPart, systemParameters);

    solveEquation(systemParameters.n * systemParameters.n);
}
