#ifndef QUARTIC_H
#define QUARTIC_H

#ifdef __cplusplus
extern "C" {
#endif
    
double real_depressed_cubic_root(double A, double B);

double real_cubic_root(double B, double C, double D);

int real_cubic_roots(double B, double C, double D,
                     double roots[]);

int real_quartic_roots(double B, double C, double D, double E,
                       double roots[]);

#ifdef __cplusplus
}
#endif
    
#endif
