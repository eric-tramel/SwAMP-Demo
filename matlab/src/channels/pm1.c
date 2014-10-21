#include "../swgamp.h"

void channel_pm1( size_t m, double *y, double *w, double *v, double *prmts, 
        double *g, double *dg, int learn ) {
    double v_eff;
    double arg, z, erfcx;
    unsigned int mu;

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
    } else {
        for (mu = 0; mu < m; mu++) channel_pm1(1, &y[mu], &w[mu], &v[mu], prmts, &g[mu], &dg[mu], 0);
    }
}
