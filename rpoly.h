#pragma once

#include <mathbits/cpp.h>

BEGIN_C_DECLS;

// Opaque state handle
struct RPoly_State;

// Allocate state for finding roots of a real polynomial.
// Degrees up to max_degree are supported by returned state.
struct RPoly_State *real_poly_alloc(int max_degree);

// Release state
void real_poly_release(struct RPoly_State *s);


// Find roots of a polynomial with real coefficients.
// p[] holds coefficients of the polynomial:
//   p(x) = p[0]*x^degree + ... + p[degree]
//
// Caller must have already allocated 'state'.
// Returns number of roots stored in zeror[], zeroi[].
int real_poly_roots_compute(const double p[], int degree,
                            struct RPoly_State *state,
                            double zeror[], double zeroi[]);

// Convenience function.
// Find roots of a polynomial with real coefficients.
// p[] holds coefficients of the polynomial:
//   p(x) = p[0]*x^degree + ... + p[degree]
//
// Returns number of roots stored in zeror[], zeroi[]
static inline int real_poly_roots(const double p[], int degree,
                                  double zeror[], double zeroi[])
{
    struct RPoly_State *state = real_poly_alloc(degree);
    int nr = real_poly_roots_compute(p, degree, state, zeror, zeroi);
    real_poly_release(state);

    return nr;
}

END_C_DECLS;
