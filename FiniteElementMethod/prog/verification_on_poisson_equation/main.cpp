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
    
    double errorSum = -10;
    double errorMax = -10;
    if (method != programMethods.RECTANGLE_SECOND_ORDER)
    {
        solveEquation(systemParameters.n * systemParameters.n);
        errorSum = compareWithAnalyticNormSum(coordinateMesh,
                                              systemParameters.n,
                                              systemParameters);
        errorMax = compareWithAnalyticNormMax(coordinateMesh,
                                              systemParameters.n,
                                              systemParameters);
    }
    else
    {
        solveEquation(systemParameters.n * systemParameters.n - ((systemParameters.n - 1) / 2) * ((systemParameters.n - 1) / 2));
        errorSum = compareWithAnalyticRectangleSecondOrderNormSum(coordinateMesh,
                                                                  systemParameters.n,
                                                                  systemParameters);
        errorMax = compareWithAnalyticRectangleSecondOrderNormMax(coordinateMesh,
                                                                  systemParameters.n,
                                                                  systemParameters);
    }

    std::cout << "Error(sum) on analytic solution of Poisson equation is: " << errorSum << std::endl;
    std::cout << "Error(max) on analytic solution of Poisson equation is: " << errorMax << std::endl;
}
