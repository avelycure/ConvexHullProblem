#include "main.hpp"

int main()
{
    std::string method;
    double *rightPart;
    Point **coordinateMesh;
    double **matrixPressure;
    SystemParameters systemParameters;

    TriangleContributionMatrix *contributionMatrix;
    TriangleRightPart *localRigthParts;

    RectangleContributionMatrix *contributionMatrixR;
    RectangleRightPart *localRigthPartsR;

    TriangleContributionMatrixSecondOrder *contributionMatrixQT;
    TriangleRightPartSecondOrder *localRigthPartsQT;

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

    if (method == "qtrig")
        solveWithTrianglesSecondOrder(contributionMatrixQT, localRigthPartsQT, coordinateMesh,
                                         matrixPressure, rightPart, systemParameters);


    solveEquation(systemParameters.n * systemParameters.n);
}
