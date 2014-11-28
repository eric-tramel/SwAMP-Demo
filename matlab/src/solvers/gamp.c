#include "../swgamp.h"

/* (Sequential) AMP for sparse matrices */
void gamp( 
        size_t n, size_t m, double *y, double *F, int *ir, int *jc, 
        void (*channel) (size_t, double*, double*, double*, double*, double*, double*, int), double *ch_prmts, int learn_channel, 
        void (*prior) (int, double*, double*, double*, double*, double*, double*, int), double *pr_prmts, int learn_prior, 
        int t_max, double eps, double damp, int disp, FILE *output, FILE *history, double *x,
        double *a, double *c, double *r, double *sig,
        int mean_removal, int calc_vfe, int adaptive_damp, int no_violations, int site_rejection  
    ) {
    double *a_proj, *c_proj, *w, *v, *g, *dg;
    double a_old, c_old, g_proj, dg_proj;
    double diff, res, mse;
    double x_n, a_n;
    int n_noaux = n-2;
    int m_noaux = m-2;

    int m_aux_a = m-2;
    int m_aux_b = m-1;
    int n_aux_a = n-2;
    int n_aux_b = n-1;

    unsigned int i, mu, idx, t;
    int *seq, key;

    /* Alloc. structures */
    a_proj = malloc(sizeof(double) * m);
    c_proj = malloc(sizeof(double) * m);
    w = malloc(sizeof(double) * m);
    v = malloc(sizeof(double) * m);
    g = malloc(sizeof(double) * m);
    dg = malloc(sizeof(double) * m);
    seq = malloc(sizeof(int) * n);

    if (!a_proj || !c_proj || !w || !v || !g || !dg || !seq)
        mexErrMsgTxt("Failure in allocating memory.");

    /* Init. variables */
    for (mu = 0; mu < m; mu++) g[mu] = 0;

    x_n = 0;
    if (x) {
        for (i = 0; i < n; i++) 
            x_n += pow(x[i], 2);
        x_n = sqrt(x_n);
    }

    /* Iterate GAMP */
    if (output) fprintf(output, "# t; MSE; RSS; diff\n");
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
                c_proj[ ir[idx] ] += (F[idx] * F[idx]) * c[i];
            }
        }

        /* Update w and v */
        for (mu = 0; mu < m; mu++) {
            w[mu] = a_proj[mu] - c_proj[mu] * g[mu];
            v[mu] = c_proj[mu];
        }
        // For Mean Removal: We want to perform the channel operation to calculate the
        //                   proper channel-based coefficients on the non-auxilliary 
        //                   terms, but we want to utilize the "impulse" channel for the
        //                   auxilliary terms themselves.
        channel((mean_removal ? m_noaux : m), y, w, v, ch_prmts, g, dg, 0);     /* Only on non-aux terms */
        if(mean_removal){
            // Double check to see if these should actually be negatives
            g[m_aux_a] = -w[m_aux_a] / v[m_aux_a];
            g[m_aux_b] = -w[m_aux_b] / v[m_aux_b];
            // Double check to see if these should actually be negatives
            dg[m_aux_a] = -1/v[m_aux_a];
            dg[m_aux_b] = -1/v[m_aux_b];
        }


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
            if(!mean_removal || (mean_removal && i < n_noaux)){
                prior(1, &r[i], &sig[i], pr_prmts, &a[i], &c[i], NULL, 0);
            }else{
                // Operate in a "priorless" mode for the auxilliary terms
                // when using mean-removal
                a[i] = r[i];
                c[i] = sig[i];
            }



            /* Update {w, v}, {g, dg} */
            for (idx = jc[i]; idx < jc[i + 1]; idx++) {
                mu = ir[idx];

                w[mu] += F[idx] * (a[i] - a_old);
                v[mu] += (F[idx] * F[idx]) * (c[i] - c_old);

                if(!mean_removal || (mean_removal && mu < m_noaux)){
                    channel(1, &y[mu], &w[mu], &v[mu], ch_prmts, &g[mu], &dg[mu], 0);
                }else{
                    // Make sure that we use the impulse prior on these auxilliary terms
                    g[mu]  = -w[mu]/v[mu];
                    dg[mu] = -1/v[mu];
                }
            }

            diff += fabs(a[i] - a_old);
        }

        /* Update channel and prior parameters */
        // Only update these terms based upon the non-auxilliary terms
        if (learn_channel) channel((mean_removal ? m_noaux : m), y, w, v, ch_prmts, g, dg, 1);
        if (learn_prior) prior((mean_removal ? n_noaux : n), r, sig, pr_prmts, a, c, NULL, 1);

        /* Compute averages */
        res = 0;
        for (mu = 0; mu < m; mu++)
            res += pow(y[mu] - a_proj[mu], 2);
        res /= m;

        mse = 0;
        if (x) {
            a_n = 0;
            // Only claculate the MSE on the non-auxilliary terms
            for (i = 0; i < (mean_removal ? n_noaux : n); i++)
                a_n += pow(a[i], 2);
            a_n = sqrt(a_n);

            for (i = 0; i < (mean_removal ? n_noaux : n); i++)
                mse += pow(a[i] / a_n - x[i] / x_n, 2);
        }

        /* Print some info. */
        if (disp)
            printf("t = %3d: mse = %.4e, RSS = %.4e, diff. = %.4e\n", t, mse, res, diff / n);
        if (output)
            fprintf(output, "%d;%g;%g;%g\n", t, mse, res, diff / n);
        mexEvalString("drawnow");

        if (history) {
            for (i = 0; i < n; i++) 
                fprintf(history, "%g ", a[i]);
            fprintf(history, "\n");
        }

        /* Check for convergence */
        if (diff / n < eps) break;
    }

    /* Dealloc. structures */
    free(seq);
    free(dg);
    free(g);
    free(v);
    free(w);
    free(c_proj);
    free(a_proj);
}
