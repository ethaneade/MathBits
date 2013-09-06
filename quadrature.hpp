#pragma once

#include <cmath>

namespace mathbits
{
    template <class Func>
    static double quad_g7k15(const Func &f, double a, double b, double &err)
    {
        static const double nw[8][3] = {
            {0.000000000000000, 0.417959183673469, 0.209482141084728},
            {0.405845151377397, 0.381830050505119, 0.190350578064785},
            {0.741531185599394, 0.279705391489277, 0.140653259715525},
            {0.949107912342759, 0.129484966168870, 0.063092092629979},
            {0.207784955007898, 0, 0.204432940075298},
            {0.586087235467691, 0, 0.169004726639267},
            {0.864864423359769, 0, 0.104790010322250},
            {0.991455371120813, 0, 0.022935322010529}
        };

        const double x0 = (a + b)*0.5;
        const double m = b - x0;
        
        double y0 = (double)f(x0);    
        double g7 = nw[0][1] * y0;
        double k15 = nw[0][2] * y0;
        
        for (int i=1; i<8; ++i) {
            double mx = m * nw[i][0];
            double yi = (double)f(x0 + mx) + (double)f(x0 - mx);
            k15 += nw[i][2] * yi;
            g7 += nw[i][1] * yi;
        }
        
        g7 *= m;
        k15 *= m;
        
        err = 200*fabs(g7 - k15);
        err *= sqrt(err);
        
        return k15;
    }


    template <class Func>
    bool integrate(const Func &f, double a, double b,
                   double &intf_a_to_b,
                   double fac_tol=1e-9, double abs_tol=1e-20)
    {
        const int max_intervals = 2000;
        double (*ab)[2] = new double[max_intervals][2];
        int n = 0;
        
        ab[0][0] = a;
        ab[0][1] = b;
        n = 1;
        double sum = 0;
        while (n > 0) {
            const double *abi = ab[--n];
            const double ai = abi[0];
            const double bi = abi[1];
            
            double err;
            double s = quad_g7k15(f, ai, bi, err);
            if (err < fac_tol*fabs(s) || err < abs_tol) {
                sum += s;
                continue;
            }

            if (n+1 >= max_intervals) {
                delete[] ab;
                return false;
            }
            
            double mid = (ai + bi)*0.5;
            ab[n][0] = ai;
            ab[n][1] = mid;
            ab[n+1][0] = mid;
            ab[n+1][1] = bi;
            n += 2;
        }
        
        delete[] ab;
        intf_a_to_b = sum;
        return true;
    }
    
}
