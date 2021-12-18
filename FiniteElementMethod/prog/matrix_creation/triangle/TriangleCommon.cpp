#include "FEMTriangles.hpp"

/**
 * Calculate area of the triangle using coordinates of its nodes
 * */
double countArea(Point pointI, Point pointJ, Point pointK)
{
    return fabs(0.5 * (pointJ.getX() * pointK.getY() - pointK.getX() * pointJ.getY() +
                       pointJ.getY() * pointI.getX() - pointK.getY() * pointI.getX() +
                       pointK.getX() * pointI.getY() - pointJ.getX() * pointI.getY()));
}