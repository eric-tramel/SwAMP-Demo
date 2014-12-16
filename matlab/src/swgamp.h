#include "mex.h"
#include <string.h>
#include <math.h>

/* SOLVERS */
void gamp ( 
        size_t n, size_t m, double *y, double *F, int *ir, int *jc, 
        void (*channel) (size_t, double*, double*, double*, double*, double*, double*, double*, int), double *ch_prmts, int learn_channel, 
        void (*prior) (int, double*, double*, double*, double*, double*, double*, int), double *pr_prmts, int learn_prior, 
        int t_max, double eps, double damp, int disp, FILE *output, FILE *history, double *x,
        double *a, double *c, double *r, double *sig,
        int mean_removal, int calc_vfe, int adaptive_damp, int no_violations, int site_rejection
    );

/* PRIORS */
void prior_gb( int n, double *r_vec, double *sig_vec, double *prmts,
        double *a, double *c, double *log_z, int learn );
void prior_binary( int n, double *r_vec, double *sig_vec, double *prmts, 
        double *a, double *c, double *log_z, int learn );

/* CHANNELS */
void channel_gaussian( size_t m, double *y, double *w, double *v, double *prmts, 
        double *g, double *dg, double *log_z, int learn );
void channel_probit( size_t m, double *y, double *w, double *v, double *prmts, 
        double *g, double *dg, double *log_z, int learn );

/* COMMON */
void sort_rand( int n, int *seq );
double gvfe( size_t n, size_t m, double *y, double *F, int *ir, int *jc,
        void (*channel) (size_t, double*, double*, double*, double*, double*, double*, double*, int), double *ch_prmts,
        double *a, double *c, double *logz_i, double *r, double *sig, double *w, double *v );

static inline double max( double a, double b ) { return a > b ? a : b; }
static inline double min( double a, double b ) { return a < b ? a : b; }
static inline int rand_int( int k ) { return rand() / (RAND_MAX / k + 1); }
