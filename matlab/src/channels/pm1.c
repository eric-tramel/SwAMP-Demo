#include "../swgamp.h"

void channel_pm1( size_t m, double *y, double *w, double *v, double *prmts, 
        double *g, double *dg, int learn ) {
    double arg;
    unsigned int mu;

    if (m == 1) {
        arg = (*y) * (*w) / sqrt(2 * (*v));

        *g = (*y) * exp(-arg * arg) / ( sqrt(.5 * M_PI * (*v)) * erfc(-arg) );
        *dg = -(*g) * ((*w) / (*v) + (*g));
    } else {
        for (mu = 0; mu < m; mu++) channel_pm1(1, &y[mu], &w[mu], &v[mu], prmts, &g[mu], &dg[mu], 0);
    }
}
