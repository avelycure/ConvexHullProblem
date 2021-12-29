//cpp libs
#include <cmath>
#include <string>
#include <fstream>
#include <iostream>

//solver
#include "solvers/gauss/GaussSystemSolver.hpp"

//input and output
#include "libs/single_include/nlohmann/json.hpp"
#include "common/classes/system/SystemParameters.hpp"
#include "common/classes/system/Methods.hpp"

//comparision with analytic solution
#include "comparision_with_analytic/Comparision.hpp"

//common funcs for all methods
#include "common/init/InitFuncs.hpp"
#include "common/classes/mesh/Point.hpp"
#include "matrix_creation/border_conditions/BorderConditions.hpp"

//triangles
#include "matrix_creation/triangle/first_order/FirstOrderTriangleFEFA.hpp"
#include "matrix_creation/triangle/first_order/FirstOrderTriangleFEGA.hpp"
#include "matrix_creation/triangle/second_order/SecondOrderTriangleFE.hpp"
#include "common/classes/contribution_matrix/triangle/TriangleRightPart.hpp"
#include "common/classes/contribution_matrix/triangle/TriangleRightPartSecondOrder.hpp"
#include "common/classes/contribution_matrix/triangle/TriangleContributionMatrixSecondOrder.hpp"
#include "common/classes/contribution_matrix/triangle/TriangleContributionMatrix.hpp"

//rectangles
#include "matrix_creation/rectangle/first_order/FirstOrderRectangleFE.hpp"
#include "matrix_creation/rectangle/second_order/SecondOrderRectangleFE.hpp"
#include "common/classes/contribution_matrix/rectangle/first_order/FirstOrderRectangleRightPart.hpp"
#include "common/classes/contribution_matrix/rectangle/first_order/FirstOrderRectangleContributionMatrix.hpp"
#include "common/classes/contribution_matrix/rectangle/second_order/SecondOrderRectangleRightPart.hpp"
#include "common/classes/contribution_matrix/rectangle/second_order/SecondOrderRectangleContributionMatrix.hpp"