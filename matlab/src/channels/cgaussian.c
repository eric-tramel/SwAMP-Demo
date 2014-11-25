#include "../cswgamp.h"

void channel_cgaussian( size_t m, complex double *y, complex double *w, double *v, double *prmts, 
        complex double *g, double *dg, int learn ) {
    double delta_n, delta_d, delta;
    unsigned int mu;

    if (learn) {
        delta_n = delta_d = 0; /* Sums: (w / v)^2 and (1 / v) */
        for (mu = 0; mu < m; mu++) {
            delta_n += pow(cabs(g[mu]), 2);
            delta_d -= dg[mu];
        }
        prmts[0] *= (delta_n / delta_d);
    }
    delta = prmts[0];

    if (m == 1) {
        *g = (*y - *w) / (delta + *v);
        *dg = -1. / (delta + *v);
    } else if (m > 1) {
        for (mu = 0; mu < m; mu++) channel_cgaussian(1, &y[mu], &w[mu], &v[mu], prmts, &g[mu], &dg[mu], 0);
    }
}
