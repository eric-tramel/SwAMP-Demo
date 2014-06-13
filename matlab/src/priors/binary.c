#include "../swamp.h"

/* Binary (0/1) */
void prior_binary( int n, double *r_vec, double *sig_vec, double *prmts, 
        double *a, double *c, double *log_z, int learn ) {
    double rho = prmts[0];
    double z, r, sig;
    double a_mean;
    int i;

    if (!learn && n == 1) {
        r = *r_vec, sig = *sig_vec;

        z = rho + (1 - rho) * exp(.5 * (1 - 2 * r) / sig);
        (*a) = rho / z;
        (*c) = (*a) * (1 - (*a));
        if (log_z != NULL)
            (*log_z) = -.5 * pow((r - 1), 2) / sig + log(z);
    } else if (!learn && n > 1) {
        for(i = 0; i < n; i++)
            prior_binary (1, &r_vec[i], &sig_vec[i], prmts, &a[i], &c[i], &log_z[i], 0);
    }

    if (learn) {
        a_mean = 0;
        for (i = 0; i < n; i++) a_mean += a[i];
        a_mean /= n;

        prmts[0] = a_mean;
    }
}
