#include "../swgamp.h"

double gvfe( size_t n, size_t m, double *y, double *F, int *ir, int *jc,
        void (*channel) (size_t, double*, double*, double*, double*, double*, double*, double*, int), double *ch_prmts,
        double *a, double *c, double *logz_i, double *r, double *sig, double *w, double *v ) {
    double *a_proj, *c_proj, *w_fp, *w_old, *v_fp, *g, *dg, *logz_mu;
    double mu_sum, i_sum;
    int i, mu, idx;

    int iter;
    double step, scale, diff;
    int iter_max = 50;
    double tol = 1e-9, reg = 1e-13;
    double test, delta;

    /* Allocate structures */
    a_proj = malloc(sizeof(double) * m);
    c_proj = malloc(sizeof(double) * m);
    w_fp = malloc(sizeof(double) * m);
    w_old = malloc(sizeof(double) * m);
    g = malloc(sizeof(double) * m);
    dg = malloc(sizeof(double) * m);
    logz_mu = malloc(sizeof(double) * m);

    if (!a_proj || !c_proj || !w_fp || !w_old || !g || !dg || !logz_mu) exit(1);

    /* Compute a_proj and c_proj */
    for (mu = 0; mu < m; mu++) a_proj[mu] = 0., c_proj[mu] = 0.;
    for (i = 0; i < n; i++) {
        for (idx = jc[i]; idx < jc[i + 1]; idx++) {
            a_proj[ ir[idx] ] += F[idx] * a[i];
            c_proj[ ir[idx] ] += (F[idx] * F[idx]) * c[i];
        }
    }

    /* Obtain fixed point \omega using Newton's method */
    for (mu = 0; mu < m; mu++) 
        w_fp[mu] = w[mu];
    v_fp = c_proj;
    for (iter = 0; iter < iter_max; iter++) {
        channel(m, y, w_fp, v_fp, ch_prmts, g, dg, NULL, 0);

        /* Newton step */
        for (mu = 0; mu < m; mu++) {
            step = w_fp[mu] - a_proj[mu] + v_fp[mu] * g[mu];
            scale = 1. + v_fp[mu] * dg[mu];
            w_old[mu] = w_fp[mu];
            w_fp[mu] -= step / (scale + reg);
        }

        /* Check for convergence */
        diff = 0.;
        for (mu = 0; mu < m; mu++) 
            diff += fabs(1. - w_fp[mu] / w_old[mu]);
        if (diff / m < tol) {
            /*printf("Finished in %d iterations. ", iter);*/
            break;
        }
    }
    /*printf ("Diff.: %.4g\n", diff / m);*/

    /*test = 0.;*/
    /*for (mu = 0; mu < m; mu++) {*/
        /*delta = ch_prmts[0];*/
        /*test += pow(w_fp[mu] - y[mu] + (v_fp[mu] + delta) * (y[mu] - a_proj[mu]) / delta, 2);*/
    /*}*/
    /*printf("%g\n", test / m);*/
    
    /* Calculate free energy */
    channel(m, y, w_fp, v_fp, ch_prmts, g, dg, logz_mu, 0);

    mu_sum = 0., i_sum = 0.;
    for (mu = 0; mu < m; mu++)
        mu_sum += logz_mu[mu] + .5 * pow(w_fp[mu] - a_proj[mu], 2) / v_fp[mu];
    for (i = 0; i < n; i++)
        i_sum += logz_i[i] + .5 * (c[i] + pow(a[i] - r[i], 2)) / sig[i];

    return -mu_sum - i_sum;
}
