#include "../swamp.h"

/* Gauss-Bernoulli */
void prior_gb( int n, double *r_vec, double *sig_vec, double *prmts,
        double *a, double *c, double *log_z, int learn ) {
    double rho = prmts[0], pr_mean = prmts[1], pr_var = prmts[2];
    double r, sig;
    
    double isv, rsc, eff, vrp, gamma;
    int i;

    if (!learn && n == 1) {
        r = *r_vec, sig = *sig_vec;

        isv = 1. / (pr_var + sig);
        rsc = .5 * (pr_mean - r) * (pr_mean - r) * isv;
        eff = (pr_mean * sig + r * pr_var) * isv;
        vrp = pr_var * sig * isv;

        /* z = rho * sqrt(vrp / pr_var) * exp(-rsc) * (1 + gamma) */
        gamma = ((1. - rho) / rho) * sqrt(pr_var / vrp) *
            exp(-.5 * r * r / sig + rsc);

        (*a) = eff / (1 + gamma);
        (*c) = max( gamma * (*a) * (*a) + vrp / (1 + gamma), 1e-19 );
        if (log_z != NULL) 
            (*log_z) = log(rho) + .5 * log(vrp / pr_var) - rsc + log1p(gamma);
    } else if (!learn && n > 1) {
        for (i = 0; i < n; i++) 
            prior_gb (1, &r_vec[i], &sig_vec[i], prmts, &a[i], &c[i], &log_z[i], 0);
    }

    if (learn) learn_prior_gb (n, r_vec, sig_vec, a, c, prmts);
}

/* Parameters learning */
void learn_prior_gb( size_t n, double *r, double *sig, 
        double *a, double *c, double *prmts ) {
    double rho = prmts[0], pr_mean = prmts[1], pr_var = prmts[2];
    double rho_n, rho_d, mean_sum, var_sum;
    double sv, cs, ex;
    int i;

    /* Sparsity */
    rho_n = rho_d = 0;
    for (i = 0; i < n; i++) {
        cs = r[i] / sig[i] + pr_mean / pr_var;
        sv = 1. / sig[i] + 1. / pr_var;
        ex = exp(.5 * cs * cs / sv - .5 * pr_mean * pr_mean / pr_var);

        rho_n += sv * a[i] / cs;
        rho_d += 1. / ( (1 - rho) + rho * ex / sqrt(pr_var * sv) );
    }
    rho = rho_n / rho_d;

    /* Mean */
    mean_sum = 0;
    for (i = 0; i < n; i++) mean_sum += a[i];
    pr_mean = mean_sum / (rho * n);

    /* Variance */
    var_sum = 0;
    for (i = 0; i < n; i++) var_sum += c[i] + a[i] * a[i];
    pr_var = var_sum / (rho * n) - pr_mean * pr_mean;

    prmts[0] = rho;
    prmts[1] = pr_mean;
    prmts[2] = pr_var;
}
