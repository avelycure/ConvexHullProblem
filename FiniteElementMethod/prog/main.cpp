#include "main.hpp"

int main()
{
    std::string method;
    Methods programMethods;
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

    if (method == programMethods.TRIANGLE_CONSTANT_HEIGHT)
        solveWithFirstOrderTriangleFEConstantHeight(contributionMatrix,
                                                    coordinateMesh,
                                                    matrixPressure,
                                                    rightPart,
                                                    systemParameters);

    if (method == programMethods.TRIANGLE_FIRST_ORDER)
        solveWithFirstOrderTriangleFE(contributionMatrix,
                                      localRigthParts,
                                      coordinateMesh,
                                      matrixPressure,
                                      rightPart,
                                      systemParameters);

    if (method == programMethods.TRIANGLE_SECOND_ORDER)
        solveWithSecondOrderTriangleFE(contributionMatrixSecondOrder,
                                       localRigthPartsSecondOrder, coordinateMesh,
                                       matrixPressure,
                                       rightPart,
                                       systemParameters);

    if (method == programMethods.RECTANGLE_FIRST_ORDER)
        solveWithFirstOrderRectangleFE(contributionMatrixRectangle,
                                       localRigthPartsRectangle,
                                       coordinateMesh,
                                       matrixPressure,
                                       rightPart,
                                       systemParameters);

    if (method == programMethods.RECTANGLE_SECOND_ORDER)
        solveWithSecondOrderRectangleFE(contributionMatrixRectangleSecondOrder,
                                        localRigthPartsRectangleSecondOrder,
                                        coordinateMesh,
                                        matrixPressure,
                                        rightPart,
                                        systemParameters);

    if (method != programMethods.RECTANGLE_SECOND_ORDER)
        solveEquation(systemParameters.n * systemParameters.n);
    else
        solveEquation(systemParameters.n * systemParameters.n - ((systemParameters.n - 1) / 2) * ((systemParameters.n - 1) / 2));
}
