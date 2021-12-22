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

    SecondOrderRectangleContributionMatrix *contributionMatrixRectangleSecondOrder;
    SecondOrderRectangleRightPart *localRigthPartsRectangleSecondOrder;

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

    if (method == "rect2")
        solveWithSecondOrderRectangleFE(contributionMatrixRectangleSecondOrder, localRigthPartsRectangleSecondOrder, coordinateMesh,
                                        matrixPressure, rightPart, systemParameters);

    int numberOfRectFE = ((systemParameters.n - 1) / 2) * ((systemParameters.n - 1) / 2);
    std::cout << systemParameters.n * systemParameters.n - numberOfRectFE << std::endl;
    if (method != "rect2")
        solveEquation(systemParameters.n * systemParameters.n);
    else
        solveEquation(systemParameters.n * systemParameters.n - numberOfRectFE);
}
