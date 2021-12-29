#include "SecondOrderTriangleFE.hpp"

void setCoefficients(double &c1, double &c2, double &c3, double &c4, double &c5, double &c6, double &c7,
                     double &c8, double &c9, double &c10, double &c11, double &c12, double &c13,
                     double &c14, double &c15, double &c16, double &c18, double &c19,
                     double s1, double s2, double s3, double s4, double s5, double s6, double s7,
                     double k1, double k2,
                     double A1, double A2, double A3, double A4,
                     double zi)
{
     //clR1
     c1 = A1 * k1 * s2 + A2 * k2 * s2 + A2 * k1 * s3 + A3 * k2 * s3 + A3 * k1 * s4 + A4 * k2 * s4 + A4 * k1 * s5 +
          A1 * s1 * (k2 - zi) - A2 * s2 * zi - A3 * s3 * zi - A4 * s4 * zi;

     //blR4
     c2 = A1 * k1 * s2 + A2 * k2 * s2 + A2 * k1 * s3 + A3 * k2 * s3 + A3 * k1 * s4 + A4 * k2 * s4 + A4 * k1 * s5 +
          A1 * s1 * (k2 - zi) - A2 * s2 * zi - A3 * s3 * zi - A4 * s4 * zi;

     //flR1
     c3 = 2.0 * A1 * k2 * s2 + 2 * A1 * k1 * s3 + 2 * A2 * k2 * s3 + 2 * A2 * k1 * s4 + 2 * A3 * k2 * s4 + 2 * A3 * k1 * s5 +
          2 * A4 * k2 * s5 + 2 * A4 * k1 * s6 - 2 * A1 * s2 * zi - 2 * A2 * s3 * zi - 2 * A3 * s4 * zi - 2 * A4 * s5 * zi;

     //clR3
     c4 = 2 * A1 * k2 * s2 + 2 * A1 * k1 * s3 + 2 * A2 * k2 * s3 + 2 * A2 * k1 * s4 +
          2 * A3 * k2 * s4 + 2 * A3 * k1 * s5 + 2 * A4 * k2 * s5 + 2 * A4 * k1 * s6 -
          2 * A1 * s2 * zi - 2 * A2 * s3 * zi - 2 * A3 * s4 * zi - 2 * A4 * s5 * zi;

     //blR2
     c5 = A1 * k2 * s2 + A1 * k1 * s3 + A2 * k2 * s3 + A2 * k1 * s4 + A3 * k2 * s4 +
          A3 * k1 * s5 + A4 * k2 * s5 + A4 * k1 * s6 - A1 * s2 * zi - A2 * s3 * zi - A3 * s4 * zi - A4 * s5 * zi;

     //elR4
     c6 = A1 * k2 * s2 + A1 * k1 * s3 + A2 * k2 * s3 + A2 * k1 * s4 + A3 * k2 * s4 + A3 * k1 * s5 +
          A4 * k2 * s5 + A4 * k1 * s6 - A1 * s2 * zi - A2 * s3 * zi - A3 * s4 * zi - A4 * s5 * zi;

     //flR3
     c7 = 4 * A1 * k2 * s3 + 4 * A1 * k1 * s4 + 4 * A2 * k2 * s4 + 4 * A2 * k1 * s5 +
          4 * A3 * k2 * s5 + 4 * A3 * k1 * s6 + 4 * A4 * k2 * s6 + 4 * A4 * k1 * s7 -
          4 * A1 * s3 * zi - 4 * A2 * s4 * zi - 4 * A3 * s5 * zi - 4 * A4 * s6 * zi;

     //dlR2
     c8 = A1 * k2 * k2 * s2 + 2 * A1 * k1 * k2 * s3 + A2 * k2 * k2 * s3 + A1 * k1 * k1 * s4 + 2 * A2 * k1 * k2 * s4 +
          A3 * k2 * k2 * s4 + A2 * k1 * k1 * s5 + 2 * A3 * k1 * k2 * s5 + A4 * k2 * k2 * s5 + A3 * k1 * k1 * s6 +
          2 * A4 * k1 * k2 * s6 + A4 * k1 * k1 * s7 - A1 * s2 * zi * zi - A2 * s3 * zi * zi -
          A3 * s4 * zi * zi - A4 * s5 * zi * zi;

     //flR2
     c9 = A1 * k2 * k2 * s2 + 2 * A1 * k1 * k2 * s3 + A2 * k2 * k2 * s3 +
          A1 * k1 * k1 * s4 + 2 * A2 * k1 * k2 * s4 + A3 * k2 * k2 * s4 +
          A2 * k1 * k1 * s5 + 2 * A3 * k1 * k2 * s5 + A4 * k2 * k2 * s5 +
          A3 * k1 * k1 * s6 + 2 * A4 * k1 * k2 * s6 + A4 * k1 * k1 * s7 -
          A1 * s2 * zi * zi - A2 * s3 * zi * zi - A3 * s4 * zi * zi - A4 * s5 * zi * zi;

     //elR5
     c10 = A1 * k2 * k2 * s2 + 2 * A1 * k1 * k2 * s3 + A2 * k2 * k2 * s3 +
           A1 * k1 * k1 * s4 + 2 * A2 * k1 * k2 * s4 + A3 * k2 * k2 * s4 +
           A2 * k1 * k1 * s5 + 2 * A3 * k1 * k2 * s5 + A4 * k2 * k2 * s5 +
           A3 * k1 * k1 * s6 + 2 * A4 * k1 * k2 * s6 + A4 * k1 * k1 * s7 -
           A1 * s2 * zi * zi - A2 * s3 * zi * zi - A3 * s4 * zi * zi - A4 * s5 * zi * zi;

     //elR3
     c11 = A1 * k2 * k2 * s2 + 2 * A1 * k1 * k2 * s3 + A2 * k2 * k2 * s3 +
           A1 * k1 * k1 * s4 + 2 * A2 * k1 * k2 * s4 + A3 * k2 * k2 * s4 +
           A2 * k1 * k1 * s5 + 2 * A3 * k1 * k2 * s5 + A4 * k2 * k2 * s5 +
           A3 * k1 * k1 * s6 + 2 * A4 * k1 * k2 * s6 + A4 * k1 * k1 * s7 -
           A1 * s2 * zi * zi - A2 * s3 * zi * zi - A3 * s4 * zi * zi - A4 * s5 * zi * zi;

     //elR1
     c13 = A1 * k1 * k2 * s2 + 0.5 * A2 * k2 * k2 * s2 + 0.5 * A1 * k1 * k1 * s3 +
           A2 * k1 * k2 * s3 + 0.5 * A3 * k2 * k2 * s3 + 0.5 * A2 * k1 * k1 * s4 +
           A3 * k1 * k2 * s4 + 0.5 * A4 * k2 * k2 * s4 + 0.5 * A3 * k1 * k1 * s5 +
           A4 * k1 * k2 * s5 + 0.5 * A4 * k1 * k1 * s6 - 0.5 * A2 * s2 * zi * zi -
           0.5 * A3 * s3 * zi * zi - 0.5 * A4 * s4 * zi * zi + 0.5 * A1 * s1 * (k2 * k2 - zi * zi);

     //clR2
     c14 = A1 * k1 * k2 * s2 + 0.5 * A2 * k2 * k2 * s2 + 0.5 * A1 * k1 * k1 * s3 +
           A2 * k1 * k2 * s3 + 0.5 * A3 * k2 * k2 * s3 + 0.5 * A2 * k1 * k1 * s4 +
           A3 * k1 * k2 * s4 + 0.5 * A4 * k2 * k2 * s4 + 0.5 * A3 * k1 * k1 * s5 +
           A4 * k1 * k2 * s5 + 0.5 * A4 * k1 * k1 * s6 - 0.5 * A2 * s2 * zi * zi -
           0.5 * A3 * s3 * zi * zi - 0.5 * A4 * s4 * zi * zi + 0.5 * A1 * s1 * (k2 * k2 - zi * zi);

     //dlR4
     c15 = 2 * A1 * k1 * k2 * s2 + A2 * k2 * k2 * s2 + A1 * k1 * k1 * s3 + 2 * A2 * k1 * k2 * s3 +
           A3 * k2 * k2 * s3 + A2 * k1 * k1 * s4 + 2 * A3 * k1 * k2 * s4 +
           A4 * k2 * k2 * s4 + A3 * k1 * k1 * s5 + 2 * A4 * k1 * k2 * s5 +
           A4 * k1 * k1 * s6 - A2 * s2 * zi * zi - A3 * s3 * zi * zi -
           A4 * s4 * zi * zi + A1 * s1 * (k2 * k2 - zi * zi);

     //blR5
     c16 = 2 * A1 * k1 * k2 * s2 + A2 * k2 * k2 * s2 + A1 * k1 * k1 * s3 + 2 * A2 * k1 * k2 * s3 +
           A3 * k2 * k2 * s3 + A2 * k1 * k1 * s4 + 2 * A3 * k1 * k2 * s4 +
           A4 * k2 * k2 * s4 + A3 * k1 * k1 * s5 + 2 * A4 * k1 * k2 * s5 +
           A4 * k1 * k1 * s6 - A2 * s2 * zi * zi - A3 * s3 * zi * zi -
           A4 * s4 * zi * zi + A1 * s1 * (k2 * k2 - zi * zi);

     //elR2
     c18 = A1 * k1 * k2 * k2 * s2 + A2 * k2 * k2 * k2 * s2 / 3.0 + A1 * k2 * s3 +
           A1 * k1 * k1 * k2 * s3 + A2 * k1 * k2 * k2 * s3 + A3 * k2 * k2 * k2 * s3 / 3.0 +
           A1 * k1 * s4 + A1 * k1 * k1 * k1 * s4 / 3.0 + A2 * k2 * s4 +
           A2 * k1 * k1 * k2 * s4 + A3 * k1 * k2 * k2 * s4 + A4 * k2 * k2 * k2 * s4 / 3.0 +
           A2 * k1 * s5 + A3 * k2 * s5 +
           k1 * (A2 * k1 * k1 + 3 * A3 * k1 * k2 + 3 * A4 * k2 * k2) * s5 / 3.0 + A3 * k1 * s6 +
           A4 * k2 * s6 + k1 * k1 * (A3 * k1 + 3 * A4 * k2) * s6 / 3.0 + A4 * k1 * s7 +
           A4 * k1 * k1 * k1 * s7 / 3.0 - A1 * s3 * zi - A2 * s4 * zi - A3 * s5 * zi -
           A4 * s6 * zi - A2 * s2 * zi * zi * zi / 3.0 - A3 * s3 * zi * zi * zi / 3.0 -
           A4 * s4 * zi * zi * zi / 3.0 + A1 * s1 * (k2 * k2 * k2 - zi * zi * zi) / 3.0;

     //dlR5
     c19 = 4 * A1 * k1 * k2 * k2 * s2 + 4 * A2 * k2 * k2 * k2 * s2 / 3.0 + 4 * A1 * k1 * k1 * k2 * s3 +
           4 * A2 * k1 * k2 * k2 * s3 + 4 * A3 * k2 * k2 * k2 * s3 / 3.0 + 4 * A1 * k1 * k1 * k1 * s4 / 3.0 +
           4 * A2 * k1 * k1 * k2 * s4 + 4 * A3 * k1 * k2 * k2 * s4 + 4 * A4 * k2 * k2 * k2 * s4 / 3.0 +
           4 * k1 * (A2 * k1 * k1 + 3 * A3 * k1 * k2 + 3 * A4 * k2 * k2) * s5 / 3.0 +
           4 * k1 * k1 * (A3 * k1 + 3 * A4 * k2) * s6 / 3.0 + 4 * A4 * k1 * k1 * k1 * s7 / 3.0 -
           4 * A2 * s2 * zi * zi * zi / 3.0 - 4 * A3 * s3 * zi * zi * zi / 3.0 - 4 * A4 * s4 * zi * zi * zi / 3.0 +
           4.0 * A1 * s1 * (k2 * k2 * k2 - zi * zi * zi) / 3.0;
}

void setFormFunctionsCoefficients(double *&a,
                                  double *&b,
                                  double *&c,
                                  double *&d,
                                  double *&e,
                                  double *&f,
                                  Point pointI,
                                  Point pointJ,
                                  Point pointK,
                                  Point pointM,
                                  Point pointN)
{
     a[0] = (pow(pointI.getX(), 2.0) * pointI.getY() * pointN.getY() +
            pointJ.getX() * pointK.getX() * pointM.getY() * pointN.getY() +
            pointI.getX() * pointI.getY() * (pointJ.getX() * (pointI.getY() - pointM.getY()) - pointN.getY() * (pointJ.getX() + pointK.getX()))) /
           ((pointI.getX() - pointJ.getX()) * (pointI.getX() - pointK.getX()) * (pointI.getY() - pointM.getY()) * (pointI.getY() - pointN.getY()));
    a[1] = (pointI.getX() * (-pointJ.getX() * pointI.getY() + pointK.getX() * pointN.getY())) /
           ((pointI.getX() - pointJ.getX()) * (pointJ.getX() - pointK.getX()) * (pointI.getY() - pointN.getY()));
    a[2] = (pointI.getX() * pointJ.getX()) / ((pointI.getX() - pointK.getX()) * (pointJ.getX() - pointK.getX()));
    a[3] = (pointI.getX() * pointI.getY()) / ((pointI.getX() - pointJ.getX()) * (pointI.getY() - pointN.getY()));
    a[4] = (pointI.getY() * pointN.getY()) / ((pointI.getY() - pointM.getY()) * (-pointM.getY() + pointN.getY()));
    a[5] = (pointI.getY() * (pointJ.getX() * pointM.getY() - pointI.getX() * pointN.getY())) /
           ((pointI.getX() - pointJ.getX()) * (pointI.getY() - pointN.getY()) * (-pointM.getY() + pointN.getY()));

    b[0] = (-pointI.getY() * (pointI.getX() + pointJ.getX()) + pointN.getY() * (pointJ.getX() + pointK.getX())) /
           ((pointI.getX() - pointJ.getX()) * (pointI.getX() - pointK.getX()) * (pointI.getY() - pointN.getY()));
    b[1] = (pointI.getY() * (pointI.getX() + pointJ.getX()) - pointN.getY() * (pointI.getX() + pointK.getX())) /
           ((pointI.getX() - pointJ.getX()) * (pointJ.getX() - pointK.getX()) * (pointI.getY() - pointN.getY()));
    b[2] = (pointI.getX() + pointJ.getX()) / ((pointI.getX() - pointK.getX()) * (-pointJ.getX() + pointK.getX()));
    b[3] = pointI.getY() / ((pointI.getX() - pointJ.getX()) * (-pointI.getY() + pointN.getY()));
    b[4] = 0.0;
    b[5] = pointI.getY() / ((pointI.getX() - pointJ.getX()) * (pointI.getY() - pointN.getY()));

    c[0] = (-pointI.getX() * (pointI.getY() + pointN.getY()) + pointJ.getX() * (pointM.getY() + pointN.getY())) /
           ((pointI.getX() - pointJ.getX()) * (pointI.getY() - pointM.getY()) * (pointI.getY() - pointN.getY()));
    c[1] = pointI.getX() / ((pointI.getX() - pointJ.getX()) * (pointI.getY() - pointN.getY()));
    c[2] = 0.0;
    c[3] = pointI.getX() / ((pointI.getX() - pointJ.getX()) * (-pointI.getY() + pointN.getY()));
    c[4] = (pointI.getY() + pointN.getY()) / ((pointI.getY() - pointM.getY()) * (pointM.getY() - pointN.getY()));
    c[5] = (-pointJ.getX() * (pointI.getY() + pointM.getY()) + pointI.getX() * (pointI.getY() + pointN.getY())) /
           ((pointI.getX() - pointJ.getX()) * (pointI.getY() - pointN.getY()) * (-pointM.getY() + pointN.getY()));

    d[0] = 1.0 / ((pointI.getX() - pointJ.getX()) * (pointI.getX() - pointK.getX()));
    d[1] = 1.0 / ((-pointI.getX() + pointJ.getX()) * (pointJ.getX() - pointK.getX()));
    d[2] = 1.0 / ((pointI.getX() - pointK.getX()) * (pointJ.getX() - pointK.getX()));
    d[3] = 0.0;
    d[4] = 0.0;
    d[5] = 0.0;

    e[0] = 1.0 / ((pointI.getX() - pointJ.getX()) * (pointI.getY() - pointN.getY()));
    e[1] = 1.0 / ((-pointI.getX() + pointJ.getX()) * (pointI.getY() - pointN.getY()));
    e[2] = 0.0;
    e[3] = 1.0 / ((pointI.getX() - pointJ.getX()) * (pointI.getY() - pointN.getY()));
    e[4] = 0.0;
    e[5] = 1.0 / ((-pointI.getX() + pointJ.getX()) * (pointI.getY() - pointN.getY()));

    f[0] = 1.0 / ((pointI.getY() - pointM.getY()) * (pointI.getY() - pointN.getY()));
    f[1] = 0.0;
    f[2] = 0.0;
    f[3] = 0.0;
    f[4] = 1.0 / ((-pointI.getY() + pointM.getY()) * (pointM.getY() - pointN.getY()));
    f[5] = 1.0 / ((pointI.getY() - pointN.getY()) * (pointM.getY() - pointN.getY()));
}