#include "../cswgamp.h"

/* (Sequential) AMP for sparse matrices */
void cgamp ( 
        size_t n, size_t m, complex double *y, complex double *F, int *ir, int *jc, 
        void (*channel) (size_t, complex double*, complex double*, double*, double*, complex double*, double*, int), double *ch_prmts, int learn_channel, 
        void (*prior) (size_t, complex double*, double*, complex double*, complex double*, double*, int), complex double *pr_prmts, int learn_prior, 
        int t_max, double eps, double damp, int disp, FILE *output, FILE *history, double *x,
        complex double *a, double *c, complex double *r, double *sig, double *diff 
    ) {
    complex double *a_proj, *w, *g, *g_old;
    complex double a_old, g_proj;
    double *c_proj, *v, *dg;
    double c_old, dg_proj, v_old;
    double res;

    unsigned int i, mu, idx, t;
    int *seq, key;

    /* Alloc. structures */
    a_proj = malloc(sizeof(complex double) * m);
    w = malloc(sizeof(complex double) * m);
    g = malloc(sizeof(complex double) * m);
    g_old = malloc(sizeof(complex double) * m);
    c_proj = malloc(sizeof(double) * m);
    v = malloc(sizeof(double) * m);
    dg = malloc(sizeof(double) * m);
    seq = malloc(sizeof(int) * n);

    if (!a_proj || !c_proj || !w || !v || !g || !g_old || !dg || !seq)
        mexErrMsgTxt("Failure in allocating memory.");

    /* Init. variables */
    for (mu = 0; mu < m; mu++) g_old[mu] = g[mu] = dg[mu] = 0;

    /* Iterate GAMP */
    if (output) fprintf(output, "# t; RSS; diff\n");
    for (t = 0; t < t_max; t++) {
        /* Generate random permutation */
        for (key = 0; key < n; key++) seq[key] = key;
        sort_rand(n, seq);

        /* Update a_proj and c_proj */
        for (mu = 0; mu < m; mu++) 
            a_proj[mu] = c_proj[mu] = 0;
        for (i = 0; i < n; i++) {
            for (idx = jc[i]; idx < jc[i + 1]; idx++) {
                a_proj[ ir[idx] ] += F[idx] * a[i];
                c_proj[ ir[idx] ] += pow(cabs(F[idx]), 2) * c[i];
            }
        }

        /* Update w and v */
        for (mu = 0; mu < m; mu++) {
            w[mu] = a_proj[mu] - c_proj[mu] * g[mu];
            v[mu] = c_proj[mu];
        }
        channel(m, y, w, v, ch_prmts, g, dg, 0);

        for (mu = 0; mu < m; mu++) g_old[mu] = g[mu];

        /* Sweep over all n variables, in random order */
        *diff = 0., res = 0.;
        for (key = 0; key < n; key++) {
            i = seq[key];

            /* Update {sig, r}, {a, c} */
            a_old = a[i], c_old = c[i];

            g_proj = dg_proj = 0.; /* Dot products: F * g and -F^2 * dg */
            for (idx = jc[i]; idx < jc[i + 1]; idx++) { 
                g_proj += conj(F[idx]) * g[ ir[idx] ];
                dg_proj -= pow(cabs(F[idx]), 2) * dg[ ir[idx] ];
            }
            sig[i] = damp * sig[i] + (1 - damp) * (1. / dg_proj);
            r[i] = damp * r[i] + (1 - damp) * (a[i] + sig[i] * g_proj);

            prior(1, &r[i], &sig[i], pr_prmts, &a[i], &c[i], 0);

            /* Update {w, v}, {g, dg} */
            for (idx = jc[i]; idx < jc[i + 1]; idx++) {
                mu = ir[idx];
                v_old = v[mu];

                v[mu] += pow(cabs(F[idx]), 2) * (c[i] - c_old);
                w[mu] += F[idx] * (a[i] - a_old);
                    - (v[mu] - v_old) * g_old[mu];

                channel(1, &y[mu], &w[mu], &v[mu], ch_prmts, &g[mu], &dg[mu], 0);
            }

            *diff += cabs(a[i] - a_old);
        }
        *diff /= n;

        /* Update channel and prior parameters */
        if (learn_channel) channel(m, y, w, v, ch_prmts, g, dg, 1);
        if (learn_prior) prior(n, r, sig, pr_prmts, a, c, 1);

        /* Compute averages */
        res = 0;
        for (mu = 0; mu < m; mu++)
            res += pow(y[mu] - cabs(a_proj[mu]), 2);
        /*res /= m;*/

        /* Print some info. */
        if (disp)
            printf("t: %3d; (y - |F * a|)^2: %.4e, diff: %.4e\n", t, res, *diff);
        if (output)
            fprintf(output, "%d;%g;%g\n", t, res, *diff);
        mexEvalString("drawnow");

        if (history) {
            for (i = 0; i < n; i++) 
                fprintf(history, "%g+%gi ", creal(a[i]), cimag(a[i]));
            fprintf(history, "\n");
        }

        /* Check for convergence */
        if (*diff < eps || isnan(*diff)) break;
    }

    /* Dealloc. structures */
    free(seq);
    free(dg);
    free(g_old);
    free(g);
    free(v);
    free(w);
    free(c_proj);
    free(a_proj);
}
