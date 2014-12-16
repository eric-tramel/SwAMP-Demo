#include "../swgamp.h"

/* (Sequential) AMP for sparse matrices */
void gamp ( 
        size_t n, size_t m, double *y, double *F, int *ir, int *jc, 
        void (*channel) (size_t, double*, double*, double*, double*, double*, double*, double*, int), double *ch_prmts, int learn_channel, 
        void (*prior) (int, double*, double*, double*, double*, double*, double*, int), double *pr_prmts, int learn_prior, 
        int t_max, double eps, double damp, int disp, FILE *output, FILE *history, double *x,
        double *a, double *c, double *r, double *sig,
        int mean_removal, int calc_vfe, int adaptive_damp, int no_violations, int site_rejection  
    ) {
    double *logz_i, *a_proj, *c_proj, *w, *v, *g, *dg, *g_old;
    double a_old, c_old, g_proj, dg_proj;
    double diff, res, mse;
    double x_n, a_n;

    int n_eff, m_eff;
    double vfe, last_vfe;

    unsigned int i, mu, idx, t;
    int *seq, key;

    /* Alloc. structures */
    logz_i = malloc(sizeof(double) * n);
    a_proj = malloc(sizeof(double) * m);
    c_proj = malloc(sizeof(double) * m);
    w = malloc(sizeof(double) * m);
    v = malloc(sizeof(double) * m);
    g = malloc(sizeof(double) * m);
    dg = malloc(sizeof(double) * m);
    g_old = malloc(sizeof(double) * m);
    seq = malloc(sizeof(int) * n);

    if (!logz_i || !a_proj || !c_proj || !w || !v || !g || !dg || !g_old || !seq)
        mexErrMsgTxt("Failure in allocating memory.");

    /* Init. variables */
    x_n = 0;
    if (x) {
        for (i = 0; i < n; i++) 
            x_n += pow(x[i], 2);
        x_n = sqrt(x_n);
    }

    n_eff = n, m_eff = m;
    if (mean_removal)
        n_eff -= 2, m_eff -= 2;

    for (mu = 0; mu < m; mu++) w[mu] = y[mu], v[mu] = 1.;
    channel(m_eff, y, w, v, ch_prmts, g, dg, NULL, 0);

    /* Iterate GAMP */
    if (output) fprintf(output, "# t; MSE; RSS; diff; VFE\n");
    for (t = 0; t < t_max; t++) {
        /* Generate random permutation */
        for (key = 0; key < n; key++) seq[key] = key;
        sort_rand(n, seq);

        /* Update a_proj and c_proj */
        for (mu = 0; mu < m; mu++) a_proj[mu] = c_proj[mu] = 0;
        for (i = 0; i < n; i++) {
            for (idx = jc[i]; idx < jc[i + 1]; idx++) {
                a_proj[ ir[idx] ] += F[idx] * a[i];
                c_proj[ ir[idx] ] += (F[idx] * F[idx]) * c[i];
            }
        }

        /* Update w and v */
        for (mu = 0; mu < m; mu++) {
            w[mu] = a_proj[mu] - c_proj[mu] * g[mu];
            v[mu] = c_proj[mu];

            g_old[mu] = g[mu];
        }
        channel(m_eff, y, w, v, ch_prmts, g, dg, NULL, 0);
        for (mu = m_eff; mu < m; mu++)
            g[mu] = -w[mu] / v[mu], dg[mu] = -1.0 / v[mu];
        /* (channel is delta for auxiliary factors) */ 

        /* Sweep over all n variables, in random order */
        diff = 0., res = 0.;
        for (key = 0; key < n; key++) {
            i = seq[key];

            /* Update {sig, r}, {a, c} */
            a_old = a[i], c_old = c[i];

            g_proj = dg_proj = 0.; /* Dot products: F * g and -F^2 * dg */
            for (idx = jc[i]; idx < jc[i + 1]; idx++) { 
                g_proj += F[idx] * g[ ir[idx] ];
                dg_proj -= (F[idx] * F[idx]) * dg[ ir[idx] ];
            }
            sig[i] = damp * sig[i] + (1 - damp) * (1. / dg_proj);
            r[i] = damp * r[i] + (1 - damp) * (a[i] + sig[i] * g_proj);

            /* Calculate a and c */
            if (i < n_eff)
                prior(1, &r[i], &sig[i], pr_prmts, &a[i], &c[i], calc_vfe ? &logz_i[i] : NULL, 0); 
            else
                a[i] = r[i], c[i] = sig[i], logz_i[i] = 1.0;
            /* (prior is uniform for auxiliary variables) */

            /* Update {w, v}, {g, dg} */
            for (idx = jc[i]; idx < jc[i + 1]; idx++) {
                mu = ir[idx];

                w[mu] += F[idx] * (a[i] - a_old)
                        - (F[idx] * F[idx]) * (c[i] - c_old) * g_old[mu];
                v[mu] += (F[idx] * F[idx]) * (c[i] - c_old);

                if (mu < m_eff)
                    channel(1, &y[mu], &w[mu], &v[mu], ch_prmts, &g[mu], &dg[mu], NULL, 0);
                else
                    g[mu] = -w[mu] / v[mu], dg[mu] = -1.0 / v[mu];
                /* (channel is delta for auxiliary factors) */ 
            }

            diff += fabs(a[i] - a_old);
        } /* End sweep */

        /* Calculate the post-sweep VFE */
        if (calc_vfe) {
            last_vfe = vfe;
            vfe = gvfe(n, m, y, F, ir, jc, channel, ch_prmts, a, c, logz_i, r, sig, w, v);
        }

        /* Update channel and prior parameters */
        if (learn_channel) channel(m_eff, y, w, v, ch_prmts, g, dg, NULL, 1);
        if (learn_prior) prior(n_eff, r, sig, pr_prmts, a, c, NULL, 1);

        /* Compute averages */
        res = 0;
        for (mu = 0; mu < m; mu++)
            res += pow(y[mu] - a_proj[mu], 2);
        res /= m;

        mse = 0;
        if (x) {
            if (channel == channel_probit) {
                a_n = 0;
                for (i = 0; i < n_eff; i++)
                    a_n += pow(a[i], 2);
                a_n = sqrt(a_n);

                for (i = 0; i < n_eff; i++) 
                    mse += pow(a[i] / a_n - x[i] / x_n, 2);
            } else {
                for (i = 0; i < n_eff; i++) 
                    mse += pow(a[i] - x[i], 2);
                mse /= n_eff;
            }
        }

        /* Print some info. */
        if (disp) {
            printf("t = %3d: mse = %.4e, RSS = %.4e, diff. = %.4e", t, mse, res, diff / n_eff);
            if (calc_vfe) printf(", vfe = %.4g", vfe); 
            printf("\n");
        }
        if (output)
            fprintf(output, "%d;%g;%g;%g;%g\n", t, mse, res, diff / n_eff, vfe);
        mexEvalString("drawnow");

        if (history) {
            for (i = 0; i < n; i++) 
                fprintf(history, "%g ", a[i]);
            fprintf(history, "\n");
        }

        /* Check for convergence */
        if (diff / n_eff < eps) break;
    }

    /* Dealloc. structures */
    free(logz_i);
    free(seq);
    free(g_old);
    free(dg);
    free(g);
    free(v);
    free(w);
    free(c_proj);
    free(a_proj);
}
