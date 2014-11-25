#include "../cswgamp.h"

void channel_cpr( size_t m, complex double *y, complex double *w, double *v, double *prmts, 
        complex double *g, double *dg, int learn ) {
    double phi, R0;
    double delta_n, delta_d, delta, damp;
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
    damp = prmts[1];

    if (m == 1) {
        phi = 2 * cabs(*y) * cabs(*w) / (delta + *v);
        R0 = gsl_sf_bessel_I1_scaled(phi) / gsl_sf_bessel_I0_scaled(phi);

        *g = damp * (*g) + (1 - damp) *
            (cabs(*y) * cexp(carg(*w) * I) * R0 - *w) / (delta + *v);
        *dg = damp * (*dg) + (1 - damp) *
            (pow(cabs(*y), 2) * (1 - R0 * R0) / pow(delta + *v, 2) - 1. / (delta + *v));
    } else if (m > 1) {
        for (mu = 0; mu < m; mu++) channel_cpr(1, &y[mu], &w[mu], &v[mu], prmts, &g[mu], &dg[mu], 0);
    }
}
