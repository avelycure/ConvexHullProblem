all:
	g++ main.cpp \
	common/init/InitFuncs.cpp \
	matrix_creation/triangle/first_order/GalerkinApproach.cpp \
	matrix_creation/triangle/first_order/FunctionalApproach.cpp \
	matrix_creation/triangle/second_order/GalerkinApproachQuadratic.cpp \
	matrix_creation/triangle/common/TriangleCommon.cpp \
	matrix_creation/rectangle/first_order/RectangleLinear.cpp \
	matrix_creation/border_conditions/BorderConditionsValues.cpp \
	matrix_creation/border_conditions/BorderConditionsDerivatives.cpp \
	solvers/gauss/GaussSystemSolver.cpp \
	-o fem_reynolds.out
	./fem_reynolds.out
