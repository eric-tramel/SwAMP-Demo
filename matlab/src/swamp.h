#include "mex.h"
#include <string.h>
#include <math.h>

/* SOLVERS */
void amp ( 
        size_t n, size_t m, double *y, double *F, int *ir, int *jc, 
        double *delta, int is_array, int learn_delta, 
        void (*prior) (int, double*, double*, double*, double*, double*, double*, int), double *prior_prmts, int learn_prior, 
        int t_max, double eps, double damp, int disp, FILE *output, FILE *history, double *x,
        double *a, double *c, double *r, double *sig,
        int mean_removal, int calc_vfe, int adaptive_damp, int no_violations, int site_rejection 
    );
void amp_alt ( 
        size_t n, size_t m, double *y, double *F, int *ir, int *jc, 
        double *delta, int is_array, int learn_delta, 
        void (*prior) (int, double*, double*, double*, double*, double*, double*, int), double *prior_prmts, int learn_prior, 
        int t_max, double eps, double damp, int disp, FILE *output, FILE *history, double *x,
        double *a, double *c, double *r, double *sig,
        int mean_removal, int calc_vfe, int adaptive_damp, int no_violations, int site_rejection  
    );
void amp_dense ( 
        size_t n, size_t m, double *y, double *F,
        double *delta, int is_array, int learn_delta, 
        void (*prior) (int, double*, double*, double*, double*, double*, double*, int), double *prior_prmts, int learn_prior, 
        int t_max, double eps, double damp, int disp, FILE *output, FILE *history, double *x, 
        double *a, double *c, double *r, double *sig,
        int mean_removal, int calc_vfe, int adaptive_damp, int no_violations, int site_rejection  
    );

/* PRIORS */
void prior_gb( int n, double *r_vec, double *sig_vec, double *prmts,
        double *a, double *c, double *log_z, int learn );
void prior_binary( int n, double *r_vec, double *sig_vec, double *prmts, 
        double *a, double *c, double *log_z, int learn );

void learn_prior_gb( size_t n, double *r, double *sig, 
        double *a, double *c, double *prmts );

/* COMMON */
void sort_rand( int n, int *seq );
double awgn_vfe (size_t n, size_t m, 
            double *y,
            double *F, int *ir, int *jc,
            double *a, double *c, double *log_z,
            double *r, double *sig, 
            double *delta, int is_array,
            int local_idx);

static inline double max( double a, double b ) { return a > b ? a : b; }
static inline double min( double a, double b ) { return a < b ? a : b; }
static inline int rand_int( int k ) { return rand() / (RAND_MAX / k + 1); }
