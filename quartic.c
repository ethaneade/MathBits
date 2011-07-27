#include <math.h>
#include <mathbits/quartic.h>

double real_depressed_cubic_root(double A, double B)
{
    static const double third = 1.0/3;

    // Cardano's method
    
    double thirdA = A * third;
    double C = thirdA * thirdA * thirdA;
    double disc = B*B + 4*C;
    
    if (disc < 0) {
        // Complex case
        double root_disc = sqrt(-disc);
        double r = sqrt(-C);
        double theta = atan2(root_disc, -B);

        // Only the real part of the cube root is needed:
        // Complex conjugation commutes with cube root, so
        // the sum of cube roots of conjugates is real
        
        double tr = cbrt(r) * cos(theta * third);
        return 2 * tr;
    } else {
        // Real case
        double root_disc = sqrt(disc);
        double t = cbrt(0.5*(-B + root_disc));
        double s = cbrt(0.5*(-B - root_disc));
        return s + t;
    }    
}

double real_cubic_root(double B, double C, double D)
{
    static const double third = 1.0/3;
    double thirdB = third * B;
    
    // Depress
    double P = C - thirdB * B;
    double Q = thirdB * (2*thirdB*thirdB - C) + D;
    
    return real_depressed_cubic_root(P, Q) - thirdB;
}

int real_cubic_roots(double B, double C, double D,
                     double roots[])
{
    // Find at least one real root
    double r0 = real_cubic_root(B, C, D);
    roots[0] = r0;

    // Divide out root
    double b = (B+r0);
    double c = C + r0*b;

    // Solve resulting quadratic
    double disc = b*b - 4*c;
    if (disc < 0) {
        return 1;
    } else if (disc == 0) {
        roots[1] = -0.5 * b;
        return 2;
    }

    double root_disc = sqrt(disc);
    roots[1] = 0.5 * (-b - root_disc);
    roots[2] = 0.5 * (-b + root_disc);
    return 3;
}


int real_quartic_roots(double B, double C, double D, double E,
                       double roots[])
{
    // Ferrari's method
    
    double fourth_B = 0.25 * B;
    double sixteenth_B_sq = fourth_B * fourth_B;
    
    double alpha, beta, gamma;
    alpha = C - 6 * sixteenth_B_sq;
    beta = D + B * (2*sixteenth_B_sq - 0.5*C);
    gamma = E + fourth_B*(fourth_B*(-3*sixteenth_B_sq + C) - D);

    double alpha2 = alpha*alpha;

    // Check for special case
    if (beta == 0) {
        double disc_inner = alpha2 - 4*gamma;
        if (disc_inner < 0)
            return 0;
        if (disc_inner == 0) {
            if (alpha > 0)
                return 0;
            else if (alpha == 0) {
                roots[0] = -fourth_B;
                return 1;
            } else {                
                double ro = sqrt(-0.5*alpha);
                roots[0] = -ro - fourth_B;
                roots[1] = ro - fourth_B;
                return 2;
            }
        }
        
        double ri = sqrt(disc_inner);
        double disc1 = -0.5*(alpha - ri);
        double disc2 = -0.5*(alpha + ri);
        
        int nr = 0;
        if (disc1 == 0) {
            roots[nr++] = -fourth_B;
        } else if (disc1 > 0) {
            double root_d1 = sqrt(disc1);
            roots[nr++] = -fourth_B - root_d1;
            roots[nr++] = -fourth_B + root_d1;
        }
        
        if (disc2 == 0) {
            roots[nr++] = -fourth_B;
        } else if (disc2 > 0) {
            double root_d2 = sqrt(disc2);
            roots[nr++] = -fourth_B - root_d2;
            roots[nr++] = -fourth_B + root_d2;
        }
        return nr;
    }

    // Nested cubic
    double b, c, d;
    b = 2.5 * alpha;
    c = 2*alpha2 - gamma;
    d = 0.5*alpha*(alpha2 - gamma) - 0.125*beta*beta;

    // Find largest real root of nested cubic
    double cr[3];
    int ncr = real_cubic_roots(b, c, d, cr);
    double y = cr[0];
    int i;
    for (i=1; i<ncr; ++i)
        y = cr[i] > y ? cr[i] : y;

    double W_sq = alpha + 2*y;
    if (W_sq <= 0)
        return -1;

    double W = sqrt(W_sq);
    
    double invW = 1/W;
    double disc_a = -3*alpha - 2*y;
    double disc_b = 2*beta*invW;

    int nr=0;
    double disc0 = disc_a + disc_b;
    if (disc0 == 0) {
        roots[nr++] = -0.5*W - fourth_B;
    } else if (disc0 > 0) {
        double rd = sqrt(disc0);
        roots[nr++] = -0.5*(W + rd) - fourth_B;
        roots[nr++] = -0.5*(W - rd) - fourth_B;
    }
    
    double disc1 = disc_a - disc_b;
    if (disc1 == 0) {
        roots[nr++] = 0.5*W - fourth_B;
    } else if (disc1 > 0) {
        double rd = sqrt(disc1);
        roots[nr++] = 0.5*(W - rd) - fourth_B;
        roots[nr++] = 0.5*(W + rd) - fourth_B;
    }
    return nr;
}

#if 0

// Simple test code

double poly_eval(const double coef[], int n, double x)
{
    double y = coef[n];
    int i;
    for (i=n-1; i>=0; --i)
        y = y*x + coef[i];
    return y;
}

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

int main(int argc, char* argv[])
{
    double coef[5];
    assert(argc == 6);
    
    int i;
    for (i=1; i<=5; ++i)
        coef[i-1] = atof(argv[i]);

    double invA = 1/coef[4];
    double roots[4];
    int nr = real_quartic_roots(invA * coef[3],
                                invA * coef[2],
                                invA * coef[1],
                                invA * coef[0],
                                roots);

    for (i=0; i<nr; ++i) {
        double x = roots[i];
        double y = poly_eval(coef, 4, x);
        printf("p(%.16g) = %.16g\n", x, y);
    }

    return 0;
}
#endif
