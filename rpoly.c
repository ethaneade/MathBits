// Adapted in 2013 by Ethan Eade from NETLIB program 493:
//    rpoly.f, Jenkins and Traub, 1972
// Original comments preserved (with small modifications) below.
// All computation happens in double precision.

#include <math.h>
#include <stdlib.h>
#include <float.h>
#include "rpoly.h"

typedef enum {false, true} bool;

#ifndef RPOLY_DEBUG
#define RPOLY_DEBUG 0
#endif

// THE MAXIMUM RELATIVE REPRESENTATION ERROR 
// WHICH CAN BE DESCRIBED AS THE SMALLEST 
// POSITIVE FLOATING POINT NUMBER SUCH THAT 
// 1.D0+ETA IS GREATER THAN 1. 
static const double fp_eta = DBL_EPSILON;

// ARE AND MRE REFER TO THE UNIT ERROR IN '+' AND '*' RESPECTIVELY. 
// THEY ARE ASSUMED TO BE THE SAME AS ETA. 
static const double fp_are = DBL_EPSILON;
static const double fp_mre = DBL_EPSILON;

struct RPoly_State
{
    int max_degree;
    
    int n;
    double *p, *qp, *k, *qk, *svk, *temp;
    double sr, si, u, v;
    double a, b, c, d;
    double a1, a2, a3, a6, a7;
    double f, g, h;
    double szr, szi, lzr, lzi;
};

struct RPoly_State *real_poly_alloc(int deg)
{
    const int n = deg + 1;
    void *buf = malloc(sizeof(double) * n * 6 + sizeof(struct RPoly_State));
    if (buf == 0)
        return 0;
    
    struct RPoly_State *s = (struct RPoly_State*)buf;
    s->p = (double*)(s + 1);
    s->qp = s->p + n;
    s->k = s->qp + n;
    s->qk = s->k + n;
    s->svk = s->qk + n;
    s->temp = s->svk + n;

    s->max_degree = deg;

    return s;
}

void real_poly_release(struct RPoly_State *s)
{
    free(s);
}

static void rpoly_quad(double, double, double,
                       double *, double *, double *, double *);

static void rpoly_quadsd(int, double, double, const double *,
                         double *, double *, double *);

static void rpoly_newest(struct RPoly_State *state, int, double *, double *);

static void rpoly_nextk(struct RPoly_State *state, int type);
static void rpoly_calcsc(struct RPoly_State *state, int *type);

static int rpoly_fixed_and_variable_shift(struct RPoly_State *state, int l2);

static int rpoly_realit(struct RPoly_State *state, double *, int *);
static int rpoly_quadit(struct RPoly_State *state, double, double);

static void copy_double(double *out, const double *in, int n)
{
    int i;
    for (i=0; i<n; ++i) {
        out[i] = in[i];
    }
}

static void maybe_scale_coefs(double p[], int n)
{
    // FIND LARGEST AND SMALLEST MODULI OF COEFFICIENTS. 
    const double big = DBL_MAX;
    const double small = DBL_MIN;
    double max_c = 0.0;
    double min_c = big;
    int i;
    for (i = 0; i <= n; ++i) {
        double x = fabs(p[i]);
        if (x > max_c) {
            max_c = x;
        }
        if (x != 0.0 && x < min_c) {
            min_c = x;
        }
    }
    
    // SCALE IF THERE ARE LARGE OR VERY SMALL COEFFICIENTS 
    double sc = small / (fp_eta * min_c);
    if (sc > 1.0) {
        if (big < max_c * sc) {
            return;
        }
    } else if (max_c < 10.0) {
        return;
    }
    
    if (sc == 0.0) {
        sc = small;
    }

    // COMPUTES A SCALE FACTOR TO MULTIPLY THE 
    // COEFFICIENTS OF THE POLYNOMIAL. THE SCALING IS DONE 
    // TO AVOID OVERFLOW AND TO AVOID UNDETECTED UNDERFLOW 
    // INTERFERING WITH THE CONVERGENCE CRITERION. 
    const double base = 2.0;
    const int power = (int)(log(sc) / log(base) + 0.5);
    if (power != 0) {
        const double factor = pow(base, power);
        for (i = 0; i <= n; ++i) {
            p[i] *= factor;
        }
    }
}


// p[] has n+1 slots
// scratch[] must also have n+1 slots
// Returns lower bound on modulus of roots of p[]
static double root_modulus_lower_bound(const double p[], int n, double scratch[])
{
    double *pt = scratch;
    int i;
    for (i = 0; i <= n; ++i) {
        pt[i] = fabs(p[i]);
    }           
    pt[n] = -pt[n];
            
    // COMPUTE UPPER ESTIMATE OF BOUND 
    double x = pow(-pt[n] / pt[0], 1.0 / n);
    if (pt[n - 1] != 0.0) {
        // IF NEWTON STEP AT THE ORIGIN IS BETTER, USE IT. 
        const double xm = -pt[n] / pt[n - 1];
        if (xm < x) {
            x = xm;
        }
    }
            
    // CHOP THE INTERVAL (0,X) BRACKETING THE ZERO 
    while (1) {
        const double xm = x * 0.1;
        double ff = pt[0];
        for (i = 1; i <= n; ++i) {
            ff = ff * xm + pt[i];
        }
        if (ff <= 0.0) {
            break;
        }
        x = xm;
    }

    // DO NEWTON ITERATION UNTIL X CONVERGES TO TWO DECIMAL PLACES 
    double dx = 0;
    do {
        double ff = pt[0];
        double df = ff;
        for (i = 1; i < n; ++i) {
            ff = ff * x + pt[i];
            df = df * x + ff;
        }
        ff = ff * x + pt[n];
        
        dx = ff / df;
        x -= dx;
    } while (fabs(dx) > 0.005 * x);
    
    return x;    
}


// FINDS THE ZEROS OF A REAL POLYNOMIAL 
// OP  - DOUBLE PRECISION VECTOR OF COEFFICIENTS IN 
//       ORDER OF DECREASING POWERS. 
// DEGREE   - INTEGER DEGREE OF POLYNOMIAL. 
// ZEROR, ZEROI - OUTPUT DOUBLE PRECISION VECTORS OF 
//                REAL AND IMAGINARY PARTS OF THE 
//                ZEROS. 
// RETURNS NUMBER OF ZEROS FOUND. 
int real_poly_roots_compute(const double *op, int degree,
                            struct RPoly_State *state,
                            double *zeror, double *zeroi)
{    
    // PEEL OFF LEADING ZEROS 
    while (degree > 0 && op[0] == 0.) {
        --degree;
        ++op;
    }

    if (degree < 1) {
        return 0;
    }

    if (state->max_degree < degree) {
        return 0;
    }
        
    state->n = degree;
    
    // REMOVE THE ZEROS AT THE ORIGIN IF ANY 
    while (op[state->n] == 0.) {
        *zeror++ = 0.0;
        *zeroi++ = 0.0;
        --state->n;
    }

    // MAKE A COPY OF THE COEFFICIENTS 
    copy_double(state->p, op, state->n + 1);

    // INITIALIZE CONSTANTS FOR SHIFT ROTATION 
    double xx = 0.70710678;
    double yy = -xx;

    while (state->n > 2) {
        // START THE ALGORITHM FOR ONE ZERO 
        int i;
        const int n = state->n;

        maybe_scale_coefs(state->p, n);
                
        // COMPUTE LOWER BOUND ON MODULI OF ZEROS. 
        const double bnd = root_modulus_lower_bound(state->p, n, state->temp);
        
        // COMPUTE THE DERIVATIVE AS THE INITIAL K POLYNOMIAL 
        // AND DO 5 STEPS WITH NO SHIFT 
        {
            const double *const p = state->p;
            double *const k = state->k;
            const double inv_n = 1.0 / (double)n;
            k[0] = p[0];
            for (i = 1; i < n; ++i) {
                k[i] = (double) (n - i) * p[i] * inv_n;
            }
            const double aa = p[n];
            const double constant_term_thresh = fabs(p[n-1]) * fp_eta * 10.0;
            bool zerok = (k[n - 1] == 0.0);
            int pass;
            for (pass = 0; pass < 5; ++pass) {
                const double cc = k[n - 1];
                if (!zerok) {
                    // USE SCALED FORM OF RECURRENCE IF VALUE OF K AT 0 IS NONZERO 
                    double t = -aa / cc;
                    for (i = n-1; i >= 1; --i) {
                        k[i] = t * k[i - 1] + p[i];
                    }
                    k[0] = p[0];
                    zerok = (fabs(k[n - 1]) <= constant_term_thresh);
                } else {
                    // USE UNSCALED FORM OF RECURRENCE 
                    for (i = n-1; i >= 1; --i) {
                        k[i] = k[i - 1];
                    }
                    k[0] = 0.0;
                    zerok = (k[n - 1] == 0.0);
                }
            }
        }
        
        // SAVE K FOR RESTARTS WITH NEW SHIFTS 
        copy_double(state->temp, state->k, n);
        
        // LOOP TO SELECT THE QUADRATIC  CORRESPONDING TO EACH NEW SHIFT 
        int nz;
        for (i = 1; i <= 20; ++i) {
            // QUADRATIC CORRESPONDS TO A DOUBLE SHIFT TO A 
            // NON-REAL POINT AND ITS COMPLEX CONJUGATE. THE POINT 
            // HAS MODULUS BND AND AMPLITUDE ROTATED BY 94 DEGREES 
            // FROM THE PREVIOUS SHIFT 
            {
                const double cosr = -0.069756474;
                const double sinr =  0.99756405;
                const double xx1 = cosr * xx - sinr * yy;
                yy = sinr * xx + cosr * yy;
                xx = xx1;
            }
            state->sr = bnd * xx;
            state->si = bnd * yy;
            state->u = state->sr * -2.0;
            state->v = bnd;
            
            // SECOND STAGE CALCULATION, FIXED QUADRATIC 
            nz = rpoly_fixed_and_variable_shift(state, i * 20);
            if (nz > 0) {
                // THE SECOND STAGE CALLS ONE OF THE THIRD STAGE 
                // ITERATIONS AND RETURNS HERE IF SUCCESSFUL. 
                // DEFLATE THE POLYNOMIAL, STORE THE ZERO(S) AND 
                // RETURN TO THE MAIN ALGORITHM. 
                *zeror++ = state->szr;
                *zeroi++ = state->szi;
                if (nz == 2) {
                    *zeror++ = state->lzr;
                    *zeroi++ = state->lzi;
                }
                state->n -= nz;
                copy_double(state->p, state->qp, state->n + 1);
                break;
            } 
        
            // IF THE ITERATION IS UNSUCCESSFUL ANOTHER QUADRATIC 
            // IS CHOSEN AFTER RESTORING K 
            copy_double(state->k, state->temp, n);
        }

        if (nz == 0) {
            // RETURN WITH FAILURE IF NO CONVERGENCE WITH 20 
            // SHIFTS 
            return degree - state->n;
        }
    }
    
    // CALCULATE THE FINAL ZERO OR PAIR OF ZEROS 
    if (state->n == 2) {
        rpoly_quad(state->p[0], state->p[1], state->p[2],
              &zeror[0], &zeroi[0],
              &zeror[1], &zeroi[1]);
    } else {
        zeror[0] = -state->p[1] / state->p[0];
        zeroi[0] = 0.;
    }
    
    return degree;
}

#if RPOLY_DEBUG
int rpoly_total_stage2;
int rpoly_total_stage3;
#endif

// COMPUTES UP TO  l2  FIXED SHIFT K-POLYNOMIALS, 
// TESTING FOR CONVERGENCE IN THE LINEAR OR QUADRATIC CASE.
// INITIATES ONE OF THE VARIABLE SHIFT ITERATIONS AND
// RETURNS THE NUMBER OF ZEROS FOUND.
static int rpoly_fixed_and_variable_shift(struct RPoly_State *state, int l2)
{
    int j;
    double s, ui, vi;
    double ss, ts, tv, vv, ots, otv, tss;
    double svu, svv;
    double tvv;
    int type;
    int iflag;
    int nz;

    double betav = 0.25;
    double betas = 0.25;
    double oss = state->sr;
    double ovv = state->v;

    const int n = state->n;
    
    // EVALUATE POLYNOMIAL BY SYNTHETIC DIVISION 
    rpoly_quadsd(n + 1, state->u, state->v, state->p,
                 state->qp, &state->a, &state->b);
    rpoly_calcsc(state, &type);
    for (j = 1; j <= l2; ++j, ovv=vv, oss=ss, otv=tv, ots=ts) {
        // CALCULATE NEXT K POLYNOMIAL AND ESTIMATE V 
	rpoly_nextk(state, type);
	rpoly_calcsc(state, &type);
	rpoly_newest(state, type, &ui, &vi);
	vv = vi;
#if RPOLY_DEBUG
        ++rpoly_total_stage2;
#endif
        // ESTIMATE S 
	ss = 0.0;
	if (state->k[n - 1] != 0.) {
	    ss = -state->p[n] / state->k[n - 1];
	}
	tv = 1.0;
	ts = 1.0;
	if (j == 1 || type == 3)
            continue;

        // COMPUTE RELATIVE MEASURES OF CONVERGENCE OF S AND V 
        // SEQUENCES 
	if (vv != 0.0) {
	    tv = fabs((vv - ovv) / vv);
	}
	if (ss != 0.0) {
	    ts = fabs((ss - oss) / ss);
	}
        // IF DECREASING, MULTIPLY TWO MOST RECENT 
        // CONVERGENCE MEASURES 
	tvv = 1.0;
	if (tv < otv) {
	    tvv = tv * otv;
	}
	tss = 1.0;
	if (ts < ots) {
	    tss = ts * ots;
	}
        // COMPARE WITH CONVERGENCE CRITERIA 
	bool vpass = tvv < betav;
	bool spass = tss < betas;
        if (!spass && !vpass)
            continue;

        // AT LEAST ONE SEQUENCE HAS PASSED THE CONVERGENCE TEST. 
        // STORE VARIABLES BEFORE ITERATING 
        svu = state->u;
        svv = state->v;
        copy_double(state->svk, state->k, n);
        s = ss;
        
        bool stry = false;
        // CHOOSE ITERATION ACCORDING TO THE FASTEST CONVERGING SEQUENCE             
        if (spass && (!vpass || tss < tvv)) {
            nz = rpoly_realit(state, &s, &iflag);
            if (nz > 0) {
                return nz;
            }
            
            // LINEAR ITERATION HAS FAILED. FLAG THAT IT HAS BEEN 
            // TRIED AND DECREASE THE CONVERGENCE CRITERION 
            stry = true;
            betas *= 0.25;
            
            if (iflag != 0) {
                // IF LINEAR ITERATION SIGNALS AN ALMOST DOUBLE REAL 
                // ZERO ATTEMPT QUADRATIC ITERATION 
                ui = -(s + s);
                vi = s * s;
                vpass = true;
            } else {
                // RESTORE VARIABLES 
                state->u = svu;
                state->v = svv;
                copy_double(state->k, state->svk, n);
            }
        }

        if (vpass) {
            nz = rpoly_quadit(state, ui, vi);
            if (nz > 0) {
                return nz;
            }
            
            // QUADRATIC ITERATION HAS FAILED. 
            // DECREASE THE CONVERGENCE CRITERION. 
            betav *= 0.25;
            
            // TRY LINEAR ITERATION IF IT HAS NOT BEEN TRIED AND 
            // THE S SEQUENCE IS CONVERGING 
            if (spass && !stry) {
                copy_double(state->k, state->svk, n);
                nz = rpoly_realit(state, &s, &iflag);
                if (nz > 0) {
                    return nz;
                }
                
                // LINEAR ITERATION HAS FAILED. FLAG THAT IT HAS BEEN 
                // TRIED AND DECREASE THE CONVERGENCE CRITERION 
                stry = true;
                betas *= 0.25;
                
                if (iflag != 0) {
                    // IF LINEAR ITERATION SIGNALS AN ALMOST DOUBLE REAL 
                    // ZERO ATTEMPT QUADRATIC ITERATION 
                    ui = -(s + s);
                    vi = s * s;
                    nz = rpoly_quadit(state, ui, vi);
                    if (nz > 0) {
                        return nz;
                    }
                    betav *= 0.25;
                }
            }
            
            // RESTORE VARIABLES 
            state->u = svu;
            state->v = svv;
            copy_double(state->k, state->svk, n);
        }
        
        // RECOMPUTE QP AND SCALAR VALUES TO CONTINUE THE 
        // SECOND STAGE 
        rpoly_quadsd(n + 1, state->u, state->v, state->p, 
                     state->qp, &state->a, &state->b);
        rpoly_calcsc(state, &type);
    }
    return 0;
}

// VARIABLE-SHIFT K-POLYNOMIAL ITERATION FOR A 
// QUADRATIC FACTOR CONVERGES ONLY IF THE ZEROS ARE 
// EQUIMODULAR OR NEARLY SO. 
// UU,VV - COEFFICIENTS OF STARTING QUADRATIC 
// RETURNS NUMBER OF ZEROS FOUND 
static int rpoly_quadit(struct RPoly_State *state, double uu, double vv)
{
    int i, j;
    double omp = 0;
    int type;
    bool tried_cluster = false;
    double relstp = 0;
    const int n = state->n;
    const int nn = n + 1;
    
    state->u = uu;
    state->v = vv;
    j = 0;
    // MAIN LOOP 
    while (1) {
        rpoly_quad(1.0, state->u, state->v,
                   &state->szr, &state->szi,
                   &state->lzr, &state->lzi);

#if RPOLY_DEBUG
        ++rpoly_total_stage3;
#endif        
        // RETURN IF ROOTS OF THE QUADRATIC ARE REAL AND NOT 
        // CLOSE TO MULTIPLE OR NEARLY EQUAL AND OF OPPOSITE 
        // SIGN 
        if (fabs(fabs(state->szr) - fabs(state->lzr)) > fabs(state->lzr) * 0.01) {
            return 0;
        }
        
        // EVALUATE POLYNOMIAL BY QUADRATIC SYNTHETIC DIVISION 
        rpoly_quadsd(nn, state->u, state->v, state->p,
                     state->qp, &state->a, &state->b);
        double mp = fabs(state->a - state->szr * state->b) + fabs(state->szi * state->b);
        
        // COMPUTE A RIGOROUS BOUND ON THE ROUNDING ERROR IN EVALUTING P 
        double ee;
        {
            const double zm = sqrt(fabs(state->v));
            const double t = -state->szr * state->b;
            ee = fabs(state->qp[0]) * 2.0;
            for (i = 1; i < n; ++i) {
                ee = ee * zm + fabs(state->qp[i]);
            }
            ee = ee * zm + fabs(state->a + t);
            ee = ((fp_mre * 5.0 + fp_are * 4.0) * ee
                  - (fp_mre * 5.0 + fp_are * 2.0) * (fabs(state->a + t) + fabs(state->b) * zm)
                  + fp_are * 2.0 * fabs(t));
        }        
        // ITERATION HAS CONVERGED SUFFICIENTLY IF THE 
        // POLYNOMIAL VALUE IS LESS THAN 20 TIMES THIS BOUND 
        if (mp <= ee * 20.0) {
            return 2;
        }
        
        ++j;
        // STOP ITERATION AFTER 20 STEPS 
        if (j > 20) {
            return 0;
        }
        
        if (j >= 2
            && !tried_cluster
            && relstp <= 0.01
            && mp >= omp)
        {
            // A CLUSTER APPEARS TO BE STALLING THE CONVERGENCE. 
            // FIVE FIXED SHIFT STEPS ARE TAKEN WITH A U,V CLOSE 
            // TO THE CLUSTER 
            if (relstp < fp_eta) {
                relstp = fp_eta;
            }
            relstp = sqrt(relstp);
            state->u -= state->u * relstp;
            state->v += state->v * relstp;
            rpoly_quadsd(nn, state->u, state->v, state->p,
                         state->qp, &state->a, &state->b);
            for (i = 0; i < 5; ++i) {
                rpoly_calcsc(state, &type);
                rpoly_nextk(state, type);
            }
            tried_cluster = true;
            j = 0;
        }
        
        omp = mp;
        
        // CALCULATE NEXT K POLYNOMIAL AND NEW U AND V 
        rpoly_calcsc(state, &type);
        rpoly_nextk(state, type);
        rpoly_calcsc(state, &type);

        double ui, vi;
        rpoly_newest(state, type, &ui, &vi);
        
        // IF VI IS ZERO THE ITERATION IS NOT CONVERGING 
        if (vi == 0.0) {
            return 0;
        }
        
        relstp = fabs((vi - state->v) / vi);
        state->u = ui;
        state->v = vi;
    }
}


// VARIABLE-SHIFT H POLYNOMIAL ITERATION FOR A REAL ZERO. 
// SSS   - STARTING ITERATE 
// RETURNS NUMBER OF ZEROS FOUND 
// IFLAG - FLAG TO INDICATE A PAIR OF ZEROS NEAR REAL AXIS. 
static int rpoly_realit(struct RPoly_State *state, double *sss, int *iflag)
{
    int i, j;
    double s, t = 0;
    double omp = 0;
    const int n = state->n;

    s = *sss;
    *iflag = 0;
    const int MAX_ITER = 11;
    for (j=1; j<=MAX_ITER; ++j) {
        // EVALUATE P AT S 
        double pv = state->p[0];
        state->qp[0] = pv;
        for (i = 1; i <= n; ++i) {
            pv = pv * s + state->p[i];
            state->qp[i] = pv;
        }
        const double mp = fabs(pv);
#if RPOLY_DEBUG
        ++rpoly_total_stage3;
#endif
        
        // COMPUTE A RIGOROUS BOUND ON THE ERROR IN EVALUATING P 
        double ee;
        {
            const double ms = fabs(s);
            ee = fp_mre / (fp_are + fp_mre) * fabs(state->qp[0]);
            for (i = 1; i <= n; ++i) {
                ee = ee * ms + fabs(state->qp[i]);
            }
            ee = (fp_are + fp_mre) * ee - fp_mre * mp;
        }
        
        // ITERATION HAS CONVERGED SUFFICIENTLY IF THE 
        // POLYNOMIAL VALUE IS LESS THAN 20 TIMES THIS BOUND 
        if (mp <= ee * 20.0) {
            state->szr = s;
            state->szi = 0.;
            return 1;
        }
        
        // STOP ITERATION AFTER 10 STEPS 
        if (j == MAX_ITER) {
            return 0;
        }
        
        if (j >= 2
            && fabs(t) <= fabs(s - t) * 0.001
            && mp > omp)
        {
            // A CLUSTER OF ZEROS NEAR THE REAL AXIS HAS BEEN ENCOUNTERED. 
            // RETURN WITH IFLAG SET TO INITIATE A QUADRATIC ITERATION 
            *iflag = 1;
            *sss = s;
            // RETURN IF THE POLYNOMIAL VALUE HAS INCREASED SIGNIFICANTLY 
            return 0;
        }
        
        omp = mp;
        
        // COMPUTE T, THE NEXT POLYNOMIAL, AND THE NEW ITERATE 
        double kv = state->k[0];
        state->qk[0] = kv;
        for (i = 1; i < n; ++i) {
            kv = kv * s + state->k[i];
            state->qk[i] = kv;
        }
        
        if (fabs(kv) <= fabs(state->k[n - 1]) * 10.0 * fp_eta) {
            // USE UNSCALED FORM 
            state->k[0] = 0.;
            for (i = 1; i < n; ++i) {
                state->k[i] = state->qk[i - 1];
            }
        } else {
            // USE THE SCALED FORM OF THE RECURRENCE IF THE VALUE 
            // OF K AT S IS NONZERO 
            t = -pv / kv;
            state->k[0] = state->qp[0];
            for (i = 1; i < n; ++i) {
                state->k[i] = t * state->qk[i - 1] + state->qp[i];
            }
        }

        kv = state->k[0];
        for (i = 1; i < n; ++i) {
            kv = kv * s + state->k[i];
        }
        
        t = 0.;
        if (fabs(kv) > fabs(state->k[n - 1]) * 10.0 * fp_eta) {
            t = -pv / kv;
        }
        s += t;
    }
    return 0;
}

// THIS ROUTINE CALCULATES SCALAR QUANTITIES USED TO 
// COMPUTE THE NEXT K POLYNOMIAL AND NEW ESTIMATES OF 
// THE QUADRATIC COEFFICIENTS. 
// TYPE - INTEGER VARIABLE SET HERE INDICATING HOW THE 
// CALCULATIONS ARE NORMALIZED TO AVOID OVERFLOW 
// SYNTHETIC DIVISION OF K BY THE QUADRATIC 1,U,V 
static void rpoly_calcsc(struct RPoly_State *state, int *type)
{
    rpoly_quadsd(state->n, state->u, state->v, state->k,
                 state->qk, &state->c, &state->d);

    const double mc = fabs(state->c);
    const double md = fabs(state->d);
    
    if (mc <= fabs(state->k[state->n - 1]) * 100.0 * fp_eta
        && md <= fabs(state->k[state->n - 2]) * 100.0 * fp_eta)
    {
        // TYPE=3 INDICATES THE QUADRATIC IS ALMOST A FACTOR OF K 
        *type = 3;
    } else if (md < mc) {
        *type = 1;
        // TYPE=1 INDICATES THAT ALL FORMULAS ARE DIVIDED BY C 
        const double inv_c = 1.0 / state->c;
        const double e = state->a * inv_c;
        state->f = state->d * inv_c;
        state->g = state->u * e;
        state->h = state->v * state->b;
        state->a3 = state->a * e + (state->h * inv_c + state->g) * state->b;
        state->a1 = state->b - state->a * state->f;
        state->a7 = state->a + state->g * state->d + state->h * state->f;
    } else {
        *type = 2;
        // TYPE=2 INDICATES THAT ALL FORMULAS ARE DIVIDED BY D 
        const double inv_d = 1.0 / state->d;
        const double e = state->a * inv_d;
        state->f = state->c * inv_d;
        state->g = state->u * state->b;
        state->h = state->v * state->b;
        state->a3 = (state->a + state->g) * e + state->h * (state->b * inv_d);
        state->a1 = state->b * state->f - state->a;
        state->a7 = (state->f + state->u) * state->a + state->h;
    }
}

// COMPUTES THE NEXT K POLYNOMIALS USING SCALARS 
// COMPUTED IN CALCSC 
static void rpoly_nextk(struct RPoly_State *state, int type)
{
    const int n = state->n;
    const double *const qp = state->qp;
    const double *const qk = state->qk;
    double *const k = state->k;
    
    if (type == 3) {
        // USE UNSCALED FORM OF THE RECURRENCE IF TYPE IS 3 
        k[0] = 0.0;
        k[1] = 0.0;
        copy_double(&k[2], &qk[0], n-2);
        return;
    }
    
    int i;
    
    const double tmp = (type == 1) ? state->b : state->a;
    if (fabs(state->a1) > fabs(tmp) * fp_eta * 10.0) {
        // USE SCALED FORM OF THE RECURRENCE 
        state->a7 /= state->a1;
        state->a3 /= state->a1;
        k[0] = qp[0];
        k[1] = qp[1] - state->a7 * qp[0];
        for (i = 2; i < n; ++i) {
            k[i] = state->a3 * qk[i - 2] - state->a7 * qp[i - 1] + qp[i];
        }
    } else {
        // IF A1 IS NEARLY ZERO THEN USE A SPECIAL FORM OF THE 
        // RECURRENCE 
        k[0] = 0.0;
        k[1] = -state->a7 * qp[0];
        for (i = 2; i < n; ++i) {
            k[i] = state->a3 * qk[i - 2] - state->a7 * qp[i - 1];
        }
    }
}

// COMPUTE NEW ESTIMATES OF THE QUADRATIC COEFFICIENTS 
// USING THE SCALARS COMPUTED IN CALCSC. 
// USE FORMULAS APPROPRIATE TO SETTING OF TYPE.     
static void rpoly_newest(struct RPoly_State *state, int type, double *uu, double *vv)
{

    if (type == 3) {
        *uu = 0.0;
        *vv = 0.0;
        return;
    }
    
   double a4, a5;
    if (type == 2) {
        a4 = (state->a + state->g) * state->f + state->h;
        a5 = (state->f + state->u) * state->c + state->v * state->d;
    } else {
        a4 = state->a + state->u * state->b + state->h * state->f;
        a5 = state->c + (state->u + state->v * state->f) * state->d;
    }
    
    // EVALUATE NEW QUADRATIC COEFFICIENTS. 
    const double inv_p1 = 1.0 / state->p[state->n];
    const double b1 = -state->k[state->n - 1] * inv_p1;
    const double b2 = -(state->k[state->n - 2] + b1 * state->p[state->n - 1]) * inv_p1;
    const double c1 = state->v * b2 * state->a1;
    const double c2 = b1 * state->a7;
    const double c3 = b1 * b1 * state->a3;
    const double c4 = c1 - c2 - c3;
    
    const double tmp = a5 + b1 * a4 - c4;
    if (tmp == 0.) {
        *uu = 0.0;
        *vv = 0.0;
    } else {
        const double inv_tmp = 1.0 / tmp;
        *uu = state->u - (state->u * (c3 + c2) + state->v * (b1 * state->a1 + b2 * state->a7)) * inv_tmp;
        *vv = state->v * (c4 * inv_tmp + 1.0);
    }
}


// DIVIDES P BY THE QUADRATIC  1,U,V  PLACING THE 
// QUOTIENT IN Q AND THE REMAINDER IN A,B     
static void rpoly_quadsd(int nn, double u, double v, const double *p,
                         double *q, double *a, double *b)
{
    // Function Body 
    double bb = p[0];
    q[0] = bb;
    
    double aa = p[1] - u * bb;
    q[1] = aa;

    int i;
    for (i = 2; i < nn; ++i) {
	double c = p[i] - u * aa - v * bb;
	q[i] = c;
	bb = aa;
	aa = c;
    }

    *a = aa;
    *b = bb;
}

// CALCULATE THE ZEROS OF THE QUADRATIC A*Z**2+B1*Z+C. 
// THE QUADRATIC FORMULA, MODIFIED TO AVOID 
// OVERFLOW, IS USED TO FIND THE LARGER ZERO IF THE 
// ZEROS ARE REAL AND BOTH ZEROS ARE COMPLEX. 
// THE SMALLER REAL ZERO IS FOUND DIRECTLY FROM THE 
// PRODUCT OF THE ZEROS C/A. 
static void rpoly_quad(double a, double b1, double c, 
                       double *sr, double *si, double *lr, double *li)
{
    if (a == 0.0) {        
        if (b1 != 0.0) {
            *sr = -(c) / b1;
        } else {
            *sr = 0.0;
        }
        *lr = 0.0;
        *si = 0.0;
        *li = 0.0;
        return;
    }

    const double inv_a = 1.0 / a;
    
    if (c == 0.0) {
        *sr = 0.0;
        *lr = -(b1) * inv_a;
        *si = 0.0;
        *li = 0.0;
        return;
    }
    
    // COMPUTE DISCRIMINANT AVOIDING OVERFLOW 
    
    const double b = b1 * 0.5;    
    const double abs_b = fabs(b);
    const double abs_c = fabs(c);
    double d, e;
    if (abs_b >= abs_c) {
        const double inv_b = 1.0 / b;
        e = 1.0 - a * inv_b * (c * inv_b);
        d = sqrt(fabs(e)) * abs_b;
    } else {
        e = a;
        if (c < 0.0) {
            e = -a;
        }
        e = b * (b / abs_c) - e;
        d = sqrt(fabs(e)) * sqrt(abs_c);
    }
    
    if (e >= 0.0) {
        // REAL ZEROS 
        if (b >= 0.0) {
            d = -d;
        }
        double num = -b + d;
        *lr = num * inv_a;
        if (num != 0.0) {
            *sr = c / num;
        } else {
            *sr = 0.0;
        }
        *si = 0.0;
        *li = 0.0;
    } else {    
        // COMPLEX CONJUGATE ZEROS 
        *sr = -b * inv_a;
        *lr = *sr;
        *si = fabs(d * inv_a);
        *li = -(*si);
    }    
}

