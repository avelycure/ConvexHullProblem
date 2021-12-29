#include "SecondOrderRectangleFE.hpp"

void setCoefficients(double &c1, double &c2, double &c3, double &c4,
                     double &c5, double &c6, double &c7, double &c8,
                     double &c9, double &c10, double &c11, double &c12,
                     double &c13, double &c14, double &c15, double &c16,
                     double &c17, double &c18, double &c19, double &c20,
                     double &c21, double &c22, double &c23, double &c24,
                     double &c25, double &c26, double &c27, double &c28,
                     double &c29, double &c30, double &c31, double &c32,
                     double &c33, double &c34, double &c35, double &c36,
                     double &c37, double &c38, double &c39, double &c40, double &c41,
                     double A1, double A2, double A3, double A4,
                     double s1, double s2, double s3, double s4, double s5, double s6, double s7, double s8,
                     double z1, double z2, double z3, double z4, double z5)
{
      //clR1
      c1 = A1 * s1 * z1 + A2 * s2 * z1 + A3 * s3 * z1 + A4 * s4 * z1;

      //blR6
      c2 = A1 * s1 * z1 + A2 * s2 * z1 + A3 * s3 * z1 + A4 * s4 * z1;

      //blR2
      c3 = A1 * s2 * z1 + A2 * s3 * z1 + A3 * s4 * z1 + A4 * s5 * z1;

      //dlR6
      c4 = A1 * s2 * z1 + A2 * s3 * z1 + A3 * s4 * z1 + A4 * s5 * z1;

      //flR1
      c5 = 2.0 * A1 * s2 * z1 + 2.0 * A2 * s3 * z1 + 2.0 * A3 * s4 * z1 + 2.0 * A4 * s5 * z1;

      //clR3
      c6 = 2.0 * A1 * s2 * z1 + 2.0 * A2 * s3 * z1 + 2.0 * A3 * s4 * z1 + 2.0 * A4 * s5 * z1;

      //blR5
      c7 = A1 * s3 * z1 + A2 * s4 * z1 + A3 * s5 * z1 + A4 * s6 * z1;

      //tlR6
      c8 = A1 * s3 * z1 + A2 * s4 * z1 + A3 * s5 * z1 + A4 * s6 * z1;

      //flR3
      c9 = 4.0 * A1 * s3 * z1 + 4.0 * A2 * s4 * z1 + 4.0 * A3 * s5 * z1 + 4.0 * A4 * s6 * z1;

      //dlR1
      c10 = A1 * s1 * z2 + A2 * s2 * z2 + A3 * s3 * z2 + A4 * s4 * z2;

      //clR2
      c11 = A1 * s1 * z2 + A2 * s2 * z2 + A3 * s3 * z2 + A4 * s4 * z2;

      //elR6
      c12 = 2.0 * A1 * s1 * z2 + 2.0 * A2 * s2 * z2 + 2.0 * A3 * s3 * z2 + 2.0 * A4 * s4 * z2;

      //flR2
      c13 = 2.0 * A1 * s2 * z2 + 2.0 * A2 * s3 * z2 + 2.0 * A3 * s4 * z2 + 2.0 * A4 * s5 * z2;

      //dlR3
      c14 = 2.0 * A1 * s2 * z2 + 2.0 * A2 * s3 * z2 + 2.0 * A3 * s4 * z2 + 2.0 * A4 * s5 * z2;

      //blR4
      c15 = 2.0 * A1 * s2 * z2 + 2.0 * A2 * s3 * z2 + 2.0 * A3 * s4 * z2 + 2.0 * A4 * s5 * z2;

      //clR5
      c16 = 2.0 * A1 * s2 * z2 + 2.0 * A2 * s3 * z2 + 2.0 * A3 * s4 * z2 + 2.0 * A4 * s5 * z2;

      //glR6
      c17 = 2.0 * A1 * s2 * z2 + 2.0 * A2 * s3 * z2 + 2.0 * A3 * s4 * z2 + 2.0 * A4 * s5 * z2;

      //dlR7
      c18 = 2 * A1 * s2 * z2 + 2 * A2 * s3 * z2 + 2 * A3 * s4 * z2 + 2 * A4 * s5 * z2;

      //tlR1
      c19 = 2.0 * A1 * s2 * z2 + 2.0 * A2 * s3 * z2 + 2.0 * A3 * s4 * z2 + 2.0 * A4 * s5 * z2;

      //elR5
      c20 = 2.0 * A1 * s3 * z2 + 2.0 * A2 * s4 * z2 + 2.0 * A3 * s5 * z2 + 2.0 * A4 * s6 * z2;

      //tlR7
      c21 = 2.0 * A1 * s3 * z2 + 2.0 * A2 * s4 * z2 + 2.0 * A3 * s5 * z2 + 2.0 * A4 * s6 * z2;

      //flR5
      c22 = 4.0 * A1 * s3 * z2 + 4.0 * A2 * s4 * z2 + 4.0 * A3 * s5 * z2 + 4.0 * A4 * s6 * z2;

      //tlR3
      c23 = 4.0 * A1 * s3 * z2 + 4.0 * A2 * s4 * z2 + 4.0 * A3 * s5 * z2 + 4.0 * A4 * s6 * z2;

      //glR1
      c24 = A1 * s1 * z3 + A2 * s2 * z3 + A3 * s3 * z3 + A4 * s4 * z3;

      //clR4
      c25 = A1 * s1 * z3 + A2 * s2 * z3 + A3 * s3 * z3 + A4 * s4 * z3;

      //dlR2
      c26 = A1 * s3 * z1 + A2 * s4 * z1 + A3 * s5 * z1 + A4 * s6 * z1 +
            A1 * s1 * z3 + A2 * s2 * z3 + A3 * s3 * z3 + A4 * s4 * z3;

      //elR7
      c27 = 4.0 * A1 * s1 * z3 + 4.0 * A2 * s2 * z3 + 4.0 * A3 * s3 * z3 + 4.0 * A4 * s4 * z3;

      //flR4
      c28 = 2.0 * A1 * s2 * z3 + 2.0 * A2 * s3 * z3 + 2.0 * A3 * s4 * z3 + 2.0 * A4 * s5 * z3;

      //dlR5
      c29 = A1 * s4 * z1 + A2 * s5 * z1 + A3 * s6 * z1 + A4 * s7 * z1 +
            2.0 * A1 * s2 * z3 + 2.0 * A2 * s3 * z3 + 2.0 * A3 * s4 * z3 + 2.0 * A4 * s5 * z3;

      //tlR2
      c30 = A1 * s4 * z1 + A2 * s5 * z1 + A3 * s6 * z1 + A4 * s7 * z1 +
            2.0 * A1 * s2 * z3 + 2.0 * A2 * s3 * z3 + 2.0 * A3 * s4 * z3 + 2.0 * A4 * s5 * z3;

      //elR4
      c31 = 4.0 * A1 * s2 * z3 + 4.0 * A2 * s3 * z3 + 4.0 * A3 * s4 * z3 + 4.0 * A4 * s5 * z3;

      //glR7
      c32 = 4.0 * A1 * s2 * z3 + 4.0 * A2 * s3 * z3 + 4.0 * A3 * s4 * z3 + 4.0 * A4 * s5 * z3;

      //glR3
      c33 = 2.0 * A1 * s2 * z3 + 2.0 * A2 * s3 * z3 + 2.0 * A3 * s4 * z3 +
            2.0 * A4 * s5 * z3;

      //tlR5
      c34 = A1 * s5 * z1 + A2 * s6 * z1 + A3 * s7 * z1 + A4 * s8 * z1 +
            4.0 * A1 * s3 * z3 + 4.0 * A2 * s4 * z3 + 4.0 * A3 * s5 * z3 + 4.0 * A4 * s6 * z3;

      //glR2
      c35 = 2.0 * A1 * s3 * z2 + 2.0 * A2 * s4 * z2 + 2.0 * A3 * s5 * z2 + 2.0 * A4 * s6 * z2 +
            A1 * s1 * z4 + A2 * s2 * z4 + A3 * s3 * z4 + A4 * s4 * z4;

      //dlR4
      c36 = 2.0 * A1 * s3 * z2 + 2.0 * A2 * s4 * z2 + 2.0 * A3 * s5 * z2 + 2.0 * A4 * s6 * z2 +
            A1 * s1 * z4 + A2 * s2 * z4 + A3 * s3 * z4 + A4 * s4 * z4;

      //glR5
      c37 = 2.0 * A1 * s4 * z2 + 2.0 * A2 * s5 * z2 + 2.0 * A3 * s6 * z2 + 2.0 * A4 * s7 * z2 +
            2.0 * A1 * s2 * z4 + 2.0 * A2 * s3 * z4 + 2.0 * A3 * s4 * z4 + 2.0 * A4 * s5 * z4;

      //tlR4
      c38 = 2.0 * A1 * s4 * z2 + 2.0 * A2 * s5 * z2 + 2.0 * A3 * s6 * z2 + 2.0 * A4 * s7 * z2 +
            2.0 * A1 * s2 * z4 + 2.0 * A2 * s3 * z4 + 2.0 * A3 * s4 * z4 + 2.0 * A4 * s5 * z4;

      //glR4
      c39 = 4.0 * A1 * s3 * z3 + 4.0 * A2 * s4 * z3 + 4.0 * A3 * s5 * z3 + 4.0 * A4 * s6 * z3 +
            A1 * s1 * z5 + A2 * s2 * z5 + A3 * s3 * z5 + A4 * s4 * z5;

      //blR7
      c40 = 2.0 * A1 * s1 * z2 + 2.0 * A2 * s2 * z2 + 2.0 * A3 * s3 * z2 + 2.0 * A4 * s4 * z2;

      //elR2
      c41 = 2.0 * A1 * s2 * z2 + 2.0 * A2 * s3 * z2 + 2.0 * A3 * s4 * z2 + 2.0 * A4 * s5 * z2;
}

void setFormFunctionsCoefficients(double *&a,
                                  double *&b,
                                  double *&c,
                                  double *&d,
                                  double *&e,
                                  double *&f,
                                  double *&g,
                                  double *&t,
                                  Point pointI,
                                  Point pointM)
{
      double xm = pointM.getX();
      double ym = pointM.getY();

      double xi = pointI.getX();
      double yi = pointI.getY();

      a[0] = (xm * ym * (xi * (-3 * yi + ym) + xm * (yi + ym))) / ((xi - xm) * (xi - xm) * (yi - ym) * (yi - ym));
      a[1] = (4 * xi * xm * ym) / ((xi - xm) * (xi - xm) * (yi - ym));
      a[2] = (xi * ym * (xm * (-3 * yi + ym) + xi * (yi + ym))) / ((xi - xm) * (xi - xm) * (yi - ym) * (yi - ym));
      a[3] = -((4 * xi * yi * ym) / ((xi - xm) * (yi - ym) * (yi - ym)));
      a[4] = (xi * yi * (xm * (yi - 3 * ym) + xi * (yi + ym))) / ((xi - xm) * (xi - xm) * (yi - ym) * (yi - ym));
      a[5] = (4 * xi * xm * yi) / ((xi - xm) * (xi - xm) * (-yi + ym));
      a[6] = (xm * yi * (xi * (yi - 3 * ym) + xm * (yi + ym))) / ((xi - xm) * (xi - xm) * (yi - ym) * (yi - ym));
      a[7] = (4 * xm * yi * ym) / ((xi - xm) * (yi - ym) * (yi - ym));

      b[0] = (ym * ((3 * xi + xm) * yi - (xi + 3 * xm) * ym)) / ((xi - xm) * (xi - xm) * (yi - ym) * (yi - ym));
      b[1] = -((4 * (xi + xm) * ym) / ((xi - xm) * (xi - xm) * (yi - ym)));
      b[2] = (ym * ((xi + 3 * xm) * yi - (3 * xi + xm) * ym)) / ((xi - xm) * (xi - xm) * (yi - ym) * (yi - ym));
      b[3] = (4 * yi * ym) / ((xi - xm) * (yi - ym) * (yi - ym));
      b[4] = (yi * (-(3 * xi + xm) * yi + (xi + 3 * xm) * ym)) / ((xi - xm) * (xi - xm) * (yi - ym) * (yi - ym));
      b[5] = (4 * (xi + xm) * yi) / ((xi - xm) * (xi - xm) * (yi - ym));
      b[6] = (yi * (-(xi + 3 * xm) * yi + (3 * xi + xm) * ym)) / ((xi - xm) * (xi - xm) * (yi - ym) * (yi - ym));
      b[7] = -((4 * yi * ym) / ((xi - xm) * (yi - ym) * (yi - ym)));

      c[0] = (xm * (xi * (3 * yi + ym) - xm * (yi + 3 * ym))) / ((xi - xm) * (xi - xm) * (yi - ym) * (yi - ym));
      c[1] = (4 * xi * xm) / ((xi - xm) * (xi - xm) * (-yi + ym));
      c[2] = (xi * (xm * (3 * yi + ym) - xi * (yi + 3 * ym))) / ((xi - xm) * (xi - xm) * (yi - ym) * (yi - ym));
      c[3] = (4 * xi * (yi + ym)) / ((xi - xm) * (yi - ym) * (yi - ym));
      c[4] = (xi * (-xi * (3 * yi + ym) + xm * (yi + 3 * ym))) / ((xi - xm) * (xi - xm) * (yi - ym) * (yi - ym));
      c[5] = (4 * xi * xm) / ((xi - xm) * (xi - xm) * (yi - ym));
      c[6] = (xm * (-xm * (3 * yi + ym) + xi * (yi + 3 * ym))) / ((xi - xm) * (xi - xm) * (yi - ym) * (yi - ym));
      c[7] = -((4 * xm * (yi + ym)) / ((xi - xm) * (yi - ym) * (yi - ym)));

      d[0] = -((xm * (yi - 5 * ym) + xi * (3 * yi + ym)) / ((xi - xm) * (xi - xm) * (yi - ym) * (yi - ym)));
      d[1] = (4 * (xi + xm)) / ((xi - xm) * (xi - xm) * (yi - ym));
      d[2] = -((xi * (yi - 5 * ym) + xm * (3 * yi + ym)) / ((xi - xm) * (xi - xm) * (yi - ym) * (yi - ym)));
      d[3] = -((4 * (yi + ym)) / ((xi - xm) * (yi - ym) * (yi - ym)));
      d[4] = -((xi * (-5 * yi + ym) + xm * (yi + 3 * ym)) / ((xi - xm) * (xi - xm) * (yi - ym) * (yi - ym)));
      d[5] = (4 * (xi + xm)) / ((xi - xm) * (xi - xm) * (-yi + ym));
      d[6] = -((xm * (-5 * yi + ym) + xi * (yi + 3 * ym)) / ((xi - xm) * (xi - xm) * (yi - ym) * (yi - ym)));
      d[7] = (4 * (yi + ym)) / ((xi - xm) * (yi - ym) * (yi - ym));

      e[0] = (2 * ym) / ((xi - xm) * (xi - xm) * (-yi + ym));
      e[1] = (4 * ym) / ((xi - xm) * (xi - xm) * (yi - ym));
      e[2] = (2 * ym) / ((xi - xm) * (xi - xm) * (-yi + ym));
      e[3] = 0.0;
      e[4] = (2 * yi) / ((xi - xm) * (xi - xm) * (yi - ym));
      e[5] = (4 * yi) / ((xi - xm) * (xi - xm) * (-yi + ym));
      e[6] = (2 * yi) / ((xi - xm) * (xi - xm) * (yi - ym));
      e[7] = 0;

      f[0] = (2 * xm) / ((-xi + xm) * (yi - ym) * (yi - ym));
      f[1] = 0.0;
      f[2] = (2 * xi) / ((xi - xm) * (yi - ym) * (yi - ym));
      f[3] = -((4 * xi) / ((xi - xm) * (yi - ym) * (yi - ym)));
      f[4] = (2 * xi) / ((xi - xm) * (yi - ym) * (yi - ym));
      f[5] = 0.0;
      f[6] = (2 * xm) / ((-xi + xm) * (yi - ym) * (yi - ym));
      f[7] = (4 * xm) / ((xi - xm) * (yi - ym) * (yi - ym));

      g[0] = 2 / ((xi - xm) * (xi - xm) * (yi - ym));
      g[1] = 4 / ((xi - xm) * (xi - xm) * (-yi + ym));
      g[2] = 2 / ((xi - xm) * (xi - xm) * (yi - ym));
      g[3] = 0.0;
      g[4] = 2 / ((xi - xm) * (xi - xm) * (-yi + ym));
      g[5] = 4 / ((xi - xm) * (xi - xm) * (yi - ym));
      g[6] = 2 / ((xi - xm) * (xi - xm) * (-yi + ym));
      g[7] = 0.0;

      t[0] = 2 / ((xi - xm) * (yi - ym) * (yi - ym));
      t[1] = 0.0;
      t[2] = -(2 / ((xi - xm) * (yi - ym) * (yi - ym)));
      t[3] = 4 / ((xi - xm) * (yi - ym) * (yi - ym));
      t[4] = -(2 / ((xi - xm) * (yi - ym) * (yi - ym)));
      t[5] = 0.0;
      t[6] = 2 / ((xi - xm) * (yi - ym) * (yi - ym));
      t[7] = -(4 / ((xi - xm) * (yi - ym) * (yi - ym)));
}