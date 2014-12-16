#include "../swgamp.h"

void channel_probit( size_t m, double *y, double *w, double *v, double *prmts, 
        double *g, double *dg, double *log_z, int learn ) {
    double delta_n, delta_d, delta, v_eff;
    double arg, z, erfcx;
    unsigned int mu;

    gsl_sf_result res;
    int err;

    if (learn) {
        delta_n = delta_d = 0; /* Sums: g^2 and -dg */
        for (mu = 0; mu < m; mu++) {
            delta_n += pow(g[mu], 2);
            delta_d -= dg[mu];
        }
        prmts[0] *= (delta_n / delta_d);
        printf("delta = %g\n", prmts[0]);
    }
    delta = prmts[0];

    if (m == 1) {
        v_eff = prmts[0] + (*v); 
        arg = -(*y) * (*w) / sqrt(2 * v_eff);
        
        /* Series expansion of scaled complementary error function.
         * Source: https://stackoverflow.com/questions/8962542 */
        if (arg < 25) {
            erfcx = exp(arg * arg) * erfc(arg);
        } else {
            z = 1. / (arg * arg);
            erfcx = 0.564189583547756287 *
                (1. + z * (-0.5 + z * (0.75 + z * (-1.875 + z * (6.5625 - 29.53125 * z))))) / arg;
        } 

        *g = (*y) / ( sqrt(.5 * M_PI * v_eff) * erfcx );
        *dg = -(*g) * ((*w) / v_eff + (*g));
        if (log_z != NULL) {
            err = gsl_sf_log_erfc_e(arg, &res);
            if (err != GSL_SUCCESS) printf("Error computing logz_mu.\n");
            *log_z = -LOG_2 + (double) res.val;
        }
    } else {
        if (log_z != NULL)
            for (mu = 0; mu < m; mu++) channel_probit(1, &y[mu], &w[mu], &v[mu], prmts, &g[mu], &dg[mu], &log_z[mu], 0);
        else
            for (mu = 0; mu < m; mu++) channel_probit(1, &y[mu], &w[mu], &v[mu], prmts, &g[mu], &dg[mu], NULL, 0);
    }
}
