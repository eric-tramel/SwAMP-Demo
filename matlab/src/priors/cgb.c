#include "../cswgamp.h"

/* Complex Gauss-Bernoulli */
void prior_cgb( size_t n, complex double *r_vec, double *sig_vec, complex double *prmts,
        complex double *a, double *c, int learn ) {
    complex double pr_mean = prmts[1];
    double rho = prmts[0], pr_var = prmts[2];
    complex double r;
    double sig;
   
    complex double eff; 
    double isv, rsc, vrp, gamma;
    int i;

    if (!learn && n == 1) {
        r = *r_vec, sig = *sig_vec;

        isv = 1. / (pr_var + sig);
        rsc = .5 * pow(cabs(pr_mean - r), 2) * isv;
        eff = (pr_mean * sig + r * pr_var) * isv;
        vrp = pr_var * sig * isv;

        /* z = 1 + gamma */
        gamma = ((1. - rho) / rho) * (pr_var / vrp) *
            exp(-.5 * pow(cabs(r), 2) / sig + rsc);

        (*a) = eff / (1 + gamma);
        (*c) = vrp / (1 + gamma) + .5 * gamma * pow(cabs(*a), 2);
    } else if (!learn && n > 1) {
        for (i = 0; i < n; i++) 
            prior_cgb (1, &r_vec[i], &sig_vec[i], prmts, &a[i], &c[i], 0);
    }

        if (learn) learn_prior_cgb (n, r_vec, sig_vec, a, c, prmts);
}

/* Parameters learning */
void learn_prior_cgb( size_t n, complex double *r, double *sig, 
        complex double *a, double *c, complex double *prmts ) {
    complex double pr_mean = prmts[1];
    double rho = creal(prmts[0]), pr_var = creal(prmts[2]);
    complex double mean_sum;
    double rho_n, rho_d, var_sum;

    double sv, cs, ex;
    int i;

    /* Sparsity */
    /*rho_n = rho_d = 0;*/
    /*for (i = 0; i < n; i++) {*/
        /*cs = cabs(r[i]) / sig[i] + cabs(pr_mean) / pr_var;*/
        /*sv = 1. / sig[i] + 1. / pr_var;*/
        /*ex = exp(.5 * pow(cabs(cs), 2) / sv */
                /*- .5 * pow(cabs(pr_mean), 2) / pr_var);*/

        /*rho_n += sv * cabs(a[i]) / cs;*/
        /*rho_d += 1. / ( (1 - rho) + rho * ex / sqrt(pr_var * sv) );*/
    /*}*/
    /*rho = rho_n / rho_d;*/

    /* Mean */
    mean_sum = 0;
    for (i = 0; i < n; i++) mean_sum += a[i];
    pr_mean = mean_sum / (rho * n);

    /* Variance */
    var_sum = 0;
    for (i = 0; i < n; i++) var_sum += c[i] + pow(cabs(a[i]), 2);
    pr_var = var_sum / (rho * n) - pow(cabs(pr_mean), 2);

    prmts[0] = rho;
    prmts[1] = pr_mean;
    prmts[2] = pr_var;
}
