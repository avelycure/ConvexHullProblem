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

    FirstOrderRectangleContributionMatrix *contributionMatrixRectangle;
    FirstOrderRectangleRightPart *localRigthPartsRectangle;

    TriangleContributionMatrixSecondOrder *contributionMatrixSecondOrder;
    TriangleRightPartSecondOrder *localRigthPartsSecondOrder;

    readSystemParameters(systemParameters, method);

    if (method == H_CONST)
        solveWithHConst(contributionMatrix, coordinateMesh, matrixPressure, rightPart,
                        systemParameters);

    if (method == H_LINEAR)
        solveWithHLinear(contributionMatrix, localRigthParts, coordinateMesh,
                         matrixPressure, rightPart, systemParameters);

    if (method == "rect")
        solveWithRectangleFiniteElements(contributionMatrixRectangle, localRigthPartsRectangle, coordinateMesh,
                                         matrixPressure, rightPart, systemParameters);

    if (method == "qtrig")
        solveWithTrianglesSecondOrder(contributionMatrixSecondOrder, localRigthPartsSecondOrder, coordinateMesh,
                                         matrixPressure, rightPart, systemParameters);


    solveEquation(systemParameters.n * systemParameters.n);
}
