#include "../header.hpp"
int main()
{
    /**
     * Linear changing of height or constant height
     * */
    string method;

    /**
     * Divide surface on triangles, mesh is coordinate of points
     * */
    Point **coordinateMesh;

    /**
     * Result of the program, values of pressure in nodes
     * */
    double **matrixPressure;

    /**
     * Values of rigth parts in local vectors of right part
     * */
    RightPart *localRigthParts;

    /**
     * Parameters of the system, values of different materials and conditions on border
     * */
    SystemPatemeters systemParameters;

    /**
     * Global stiffness matrix
     * */
    ContributionMatrix *contributionMatrix;

    readSystemParameters(systemParameters, method);
    dimensionlessSystemParameters(systemParameters, method);

    if (method == H_CONST)
        solveWithHConst(contributionMatrix, coordinateMesh, matrixPressure,
                        systemParameters);

    if (method == H_LINEAR)
        solveWithHLinear(contributionMatrix, localRigthParts, coordinateMesh,
                         matrixPressure, systemParameters);

    const int MATRIX_PRESSURE_SIZE = systemParameters.n * systemParameters.n;
    double *rightPart = new double[MATRIX_PRESSURE_SIZE];

    if (systemParameters.borderConditions == "common")
        addBorderConditions(matrixPressure, systemParameters.n, MATRIX_PRESSURE_SIZE,
                            systemParameters.LOW_BORDER, systemParameters.HIGH_BORDER);

    if (systemParameters.borderConditions == "user")
    {
        double *borderValues = new double[MATRIX_PRESSURE_SIZE];
        cout << "in user branch" << endl;

        inputVector("border_conditions/right.txt", borderValues, MATRIX_PRESSURE_SIZE);
        cout << borderValues[1] << endl;
        addBorderConditions(matrixPressure, borderValues, systemParameters.n, MATRIX_PRESSURE_SIZE, "RIGHT");

        inputVector("border_conditions/top.txt", borderValues, MATRIX_PRESSURE_SIZE);
        cout << borderValues[1] << endl;
        addBorderConditions(matrixPressure, borderValues, systemParameters.n, MATRIX_PRESSURE_SIZE, "TOP");

        inputVector("border_conditions/left.txt", borderValues, MATRIX_PRESSURE_SIZE);
        cout << borderValues[1] << endl;
        addBorderConditions(matrixPressure, borderValues, systemParameters.n, MATRIX_PRESSURE_SIZE, "LEFT");

        inputVector("border_conditions/bottom.txt", borderValues, MATRIX_PRESSURE_SIZE);
        cout << borderValues[1] << endl;
        addBorderConditions(matrixPressure, borderValues, systemParameters.n, MATRIX_PRESSURE_SIZE, "BOTTOM");
    }

    if (systemParameters.borderConditions == "der")
    {
        cout << "in der branch1" << endl;
        double h = systemParameters.L / (systemParameters.n - 1);
        addBorderConditionsToLeftAndRight(matrixPressure,
                                          systemParameters.n,
                                          h,
                                          MATRIX_PRESSURE_SIZE,
                                          systemParameters.HIGH_BORDER,
                                          systemParameters.LOW_BORDER);
    }

    if (systemParameters.borderConditions == "der2")
    {
        cout << "in der branch2" << endl;
        double h = systemParameters.L / (systemParameters.n - 1);
        addBorderConditionsSecondOrder(matrixPressure,
                                          systemParameters.n,
                                          h,
                                          MATRIX_PRESSURE_SIZE,
                                          systemParameters.HIGH_BORDER,
                                          systemParameters.LOW_BORDER);
    }
    
    outputPressureMatrix(matrixPressure, MATRIX_PRESSURE_SIZE);

    solveEquation(systemParameters.n * systemParameters.n);
}