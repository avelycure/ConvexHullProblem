all:
	g++ main.cpp \
	common/init/InitFuncs.cpp \
	matrix_creation/triangle/first_order/FirstOrderTriangleFEGA.cpp \
	matrix_creation/triangle/first_order/FirstOrderTriangleFEFA.cpp \
	matrix_creation/triangle/second_order/SecondOrderTriangleFE.cpp \
	matrix_creation/triangle/second_order/SecondOrderTriangleCoefficients.cpp \
	matrix_creation/triangle/common/TriangleCommon.cpp \
	matrix_creation/rectangle/first_order/FirstOrderRectangleFE.cpp \
	matrix_creation/rectangle/second_order/SecondOrderRectangleFE.cpp \
	matrix_creation/rectangle/second_order/SecondOrderRectangleCoefficients.cpp \
	matrix_creation/border_conditions/BorderConditionsValues.cpp \
	matrix_creation/border_conditions/BorderConditionsDerivatives.cpp \
	solvers/gauss/GaussSystemSolver.cpp \
	comparision_with_analytic/Comparision.cpp \
	-o fem_poisson.out
	./fem_poisson.out
