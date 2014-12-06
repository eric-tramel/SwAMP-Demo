#include "../swgamp.h"

void channel_gaussian( size_t m, double *y, double *w, double *v, double *prmts, 
        double *g, double *dg, double *log_z, int learn ) {
    double delta_n, delta_d, delta;
    unsigned int mu;

    if (learn) {
        delta_n = delta_d = 0; /* Sums: (w / v)^2 and (1 / v) */
        for (mu = 0; mu < m; mu++) {
            delta_n += pow(g[mu], 2);
            delta_d -= dg[mu];
        }
        prmts[0] *= (delta_n / delta_d);
    }
    delta = prmts[0];

    if (m == 1) {
        *g = (*y - *w) / (delta + *v);
        *dg = -1. / (delta + *v);
        if (log_z != NULL)
            *log_z = -.5 * log(delta + *v) - .5 * pow(*y - *w, 2) / (delta + *v);
    } else if (m > 1) {
        if (log_z != NULL)
            for (mu = 0; mu < m; mu++) channel_gaussian(1, &y[mu], &w[mu], &v[mu], prmts, &g[mu], &dg[mu], &log_z[mu], 0);
        else
            for (mu = 0; mu < m; mu++) channel_gaussian(1, &y[mu], &w[mu], &v[mu], prmts, &g[mu], &dg[mu], NULL, 0);
    }
}
