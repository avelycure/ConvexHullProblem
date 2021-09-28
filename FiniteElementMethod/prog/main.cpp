#include "header.hpp"
int main()
{
    /**
     * Linear changing of height or constant height
     * */
    string method;

    /**
     * Vector of right parts
     * */
    double *rightPart;

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

    if (method == H_CONST)
        solveWithHConst(contributionMatrix, coordinateMesh, matrixPressure,
                        systemParameters);

    if (method == H_LINEAR)
        solveWithHLinear(contributionMatrix, localRigthParts, coordinateMesh,
                         matrixPressure, rightPart, systemParameters);
                         
    solveEquation(systemParameters.n * systemParameters.n);
}