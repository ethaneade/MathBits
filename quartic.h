#ifndef MATHBITS_QUARTIC_H
#define MATHBITS_QUARTIC_H

#include <mathbits/cpp.h>

BEGIN_C_DECLS;

// Returns one real root of the cubic x^3 + Ax + B
double real_depressed_cubic_root(double A, double B);

// Returns one real root of the cubic x^3 + Bx^2 + Cx + D
double real_cubic_root(double B, double C, double D);

// Find all real roots of the cubic x^3 + Bx^2 + Cx + D
// Stores up to 3 abscissas in roots[].
// Returns number of roots found.
int real_cubic_roots(double B, double C, double D,
                     double roots[]);

// Find all real roots of the quartic x^4 + Bx^3 + Cx^2 + Dx + E
// Stores up to 4 abscissas in roots[].
// Returns number of roots found, or -1 on numerical error.
int real_quartic_roots(double B, double C, double D, double E,
                       double roots[]);

END_C_DECLS;

#endif
