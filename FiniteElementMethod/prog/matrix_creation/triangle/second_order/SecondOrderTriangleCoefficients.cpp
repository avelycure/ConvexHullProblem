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

    //dlR5 ???
    c19 = 4 * A1 * k1 * k2 * k2 * s2 + 4 * A2 * k2 * k2 * k2 * s2 / 3.0 + 4 * A1 * k1 * k1 * k2 * s3 +
          4 * A2 * k1 * k2 * k2 * s3 + 4 * A3 * k2 * k2 * k2 * s3 / 3.0 + 4 * A1 * k1 * k1 * k1 * s4 / 3.0 +
          4 * A2 * k1 * k1 * k2 * s4 + 4 * A3 * k1 * k2 * k2 * s4 + 4 * A4 * k2 * k2 * k2 * s4 / 3.0 +
          4 * k1 * (A2 * k1 * k1 + 3 * A3 * k1 * k2 + 3 * A4 * k2 * k2) * s5 / 3.0 +
          4 * k1 * k1 * (A3 * k1 + 3 * A4 * k2) * s6 / 3.0 + 4 * A4 * k1 * k1 * k1 * s7 / 3.0 -
          4 * A2 * s2 * zi * zi * zi / 3.0 - 4 * A3 * s3 * zi * zi * zi / 3.0 - 4 * A4 * s4 * zi * zi * zi / 3.0 +
          4.0 * A1 * s1 * (k2 * k2 * k2 - zi * zi * zi) / 3.0;
}