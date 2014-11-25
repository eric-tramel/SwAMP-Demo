#include "mex.h"
#include <string.h>
#include <math.h>
#include <complex.h>
#include <gsl/gsl_sf_bessel.h>

/* SOLVERS */
void cgamp ( 
        size_t n, size_t m, complex double *y, complex double *F, int *ir, int *jc, 
        void (*channel) (size_t, complex double*, complex double*, double*, double*, complex double*, double*, int), double *ch_prmts, int learn_channel, 
        void (*prior) (size_t, complex double*, double*, complex double*, complex double*, double*, int), complex double *pr_prmts, int learn_prior, 
        int t_max, double eps, double damp, int disp, FILE *output, FILE *history, double *x,
        complex double *a, double *c, complex double *r, double *sig, double *diff
    );

/* PRIORS */
void prior_cgb( size_t n, complex double *r_vec, double *sig_vec, complex double *prmts,
        complex double *a, double *c, int learn );

/* CHANNELS */
void channel_cgaussian( size_t m, complex double *y, complex double *w, double *v, double *prmts, 
        complex double *g, double *dg, int learn );
void channel_cpr( size_t m, complex double *y, complex double *w, double *v, double *prmts, 
        complex double *g, double *dg, int learn );

void learn_prior_cgb( size_t n, complex double *r, double *sig, 
        complex double *a, double *c, complex double *prmts );

/* COMMON */
void sort_rand( int n, int *seq );

static inline double max( double a, double b ) { return a > b ? a : b; }
static inline double min( double a, double b ) { return a < b ? a : b; }
static inline int rand_int( int k ) { return rand() / (RAND_MAX / k + 1); }
