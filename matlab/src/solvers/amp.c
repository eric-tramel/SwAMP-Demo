#include "../swamp.h"

/* (Sequential) AMP for sparse matrices */
void amp ( 
        size_t n, size_t m, double *y, double *F, int *ir, int *jc, 
        double *delta, int is_array, int learn_delta, 
        void (*prior) (int, double*, double*, double*, double*, double*, double*, int), double *prior_prmts, int learn_prior, 
        int t_max, double eps, double damp, int disp, FILE *output, FILE *history, double *x,
        double *a, double *c, double *r, double *sig,
        int mean_removal, int calc_vfe, int adaptive_damp, int no_violations, int site_rejection  
    ) {
    double *w_r, *v, *g, *a_proj, *c_proj, *delta0, *logz;
    double w_proj, v_proj, a_old, c_old, r_old, s_old, logz_old, diff, res;
    double site_violation_count = 0;
    double site_violation = 0;
    double delta_n, delta_d, delta_mean, gamma;
    double mse;
    double vfe = 0;
    double last_vfe = 0;
    double local_vfe_fluctuation;
    double a_tmp,c_tmp,logz_tmp,r_tmp,sig_tmp;
    int n_noaux = n-2;
    int m_noaux = m-2;

    /* Hard coding adaptive damping params for now */
    double damp_max = 0.999;
    double damp_min = 0.0001;
    double damp_modifier_up   = 0.3;
    double damp_modifier_down = 0.08;
    double damp_site_violation = 0.9;

    unsigned int i, mu, idx, t;
    int *seq, key;

    /* Alloc. structures */
    if(calc_vfe){
        // Only allocating the space for logz if we are actually going to use it.
        logz = malloc(sizeof(double) * n);
    }
    w_r = malloc(sizeof(double) * m);
    v = malloc(sizeof(double) * m);
    g = malloc(sizeof(double) * m);
    a_proj = malloc(sizeof(double) * m);
    c_proj = malloc(sizeof(double) * m);
    delta0 = malloc(sizeof(double) * m);
    seq = malloc(sizeof(int) * n);
    if (!w_r || !v || !g || !a_proj || !c_proj || !delta0 || !seq)
        mexErrMsgTxt("Failure in allocating memory.");

    /* Init. variables */
    for (mu = 0; mu < m; mu++) delta0[mu] = delta[mu];
    for (mu = 0; mu < m; mu++) w_r[mu] = 0.;
    for (mu = 0; mu < m; mu++) v[mu] = 1.;


    /* Iterate AMP */
    if (output)
        fprintf(output, "'#t';'mse';'delta';'RSS';'diff';'F'\n");
    for (t = 0; t < t_max; t++) {
        /* Generate random permutation */
        for (key = 0; key < n; key++) seq[key] = key;
        sort_rand(n, seq);


        /* Update a_proj and c_proj */
        for (mu = 0; mu < m; mu++) a_proj[mu] = c_proj[mu] = 0;
        for (i = 0; i < n; i++)
            for (idx = jc[i]; idx < jc[i + 1]; idx++) {
                a_proj[ ir[idx] ] += F[idx] * a[i];
                c_proj[ ir[idx] ] += (F[idx] * F[idx]) * c[i];
            }

        /* Update w_r and v */
        for (mu = 0; mu < m; mu++) {
            g[mu] = w_r[mu] / v[mu];
            w_r[mu] = (y[mu] - a_proj[mu]) + c_proj[mu] * g[mu];
            
            if(mean_removal && mu >= (m - 2)){
                // Use a reall small noise to simulate the delta/impulse channel for the 
                // auxiliary coefficients
                v[mu] = c_proj[mu];
            } else {
                v[mu] = (is_array ? delta[mu] : *delta) + c_proj[mu];
            }
        }

        /* Sweep over all n variables, in random order */
        diff = res = 0.;
        site_violation_count = 0;
        for (key = 0; key < n; key++) {
            i = seq[key];
            a_old = a[i], c_old = c[i];

            /* Update r and sig ... */
            w_proj = v_proj = 0.; /* Dot products: Fw_r / v and F^2 / v */
            for (idx = jc[i]; idx < jc[i + 1]; idx++) { 
                w_proj += F[idx] * w_r[ ir[idx] ] / v[ ir[idx] ];
                v_proj += (F[idx] * F[idx]) / v[ ir[idx] ];
            }

            sig_tmp = damp * sig[i] + (1 - damp) * (1. / v_proj);
            r_tmp = damp * r[i] + (1 - damp) * (a[i] + sig_tmp * w_proj);

            /* ... then, a and c ... */
            if (calc_vfe) {
                prior(1, &r_tmp, &sig_tmp, prior_prmts, &a_tmp, &c_tmp, &logz_tmp, 0);                
            } else {
                prior(1, &r_tmp, &sig_tmp, prior_prmts, &a_tmp, &c_tmp, NULL, 0);
            }

            // If we are in mean-removal mode and this is an auxiliary variable, set
            // {a,c} to be equal to {r,sigma}
            if (mean_removal && i >= n_noaux) {
                a_tmp = r_tmp;
                c_tmp = sig_tmp;
            }

            // Don't check for site rejections if we are on an auxiliary variable in mean_removal mode
            if ((mean_removal && i < n_noaux) || !mean_removal) {
                /* Check for a local violation */
                if(site_rejection){
                    /* What are the local characteristics of the VFE when chaning this one variable ?*/ 
                    local_vfe_fluctuation = awgn_vfe_sitechange (n,m,y,
                                                                F,ir,jc,
                                                                a_tmp,c_tmp,logz_tmp,r_tmp,sig_tmp,
                                                                a,c,logz,r,sig,
                                                                delta,is_array,
                                                                i);
                    // printf("    [%d] Local VFE Fluctuation : %g\n",i,local_vfe_fluctuation);
                    if(local_vfe_fluctuation < 0){
                        /* Apply a strong damping to site */
                        a_tmp = damp_site_violation*a[i] + (1-damp_site_violation)*a_tmp;
                        c_tmp = damp_site_violation*c[i] + (1-damp_site_violation)*c_tmp;

                        site_violation_count += 1;
                        site_violation = 1;                    
                    }else{
                        site_violation = 0;
                    }
                }
            }

            // Undo a site violation if this is just the first sweep
            if(t < 1 && site_violation) site_violation = 0;

            // Now, update the state vectors if we are in the proper modes
            if(!site_rejection || !site_violation || !no_violations){
                /* ... and finally, w_r and v. */
                for (idx = jc[i]; idx < jc[i + 1]; idx++) {
                    w_r[ ir[idx] ] += F[idx] * (a[i] - a_tmp)
                        - (F[idx] * F[idx]) * (c[i] - c_tmp) * g[ ir[idx] ];
                    v[ ir[idx] ] += (F[idx] * F[idx]) * (c_tmp - c[i]);
                }

                diff += fabs(a_tmp - a[i]);
            
                a[i] = a_tmp;
                c[i] = c_tmp;                
                r[i] = r_tmp;
                sig[i] = sig_tmp;    

                // printf("Updated Site %d\n",i);

                if(calc_vfe) logz[i] = logz_tmp;        /* Otherwise we may attempt to write to the unallocated array */
            }
        }   /* End Sweep */

        /* Calculate the post-sweep VFE */
        if (calc_vfe){
            last_vfe = vfe;
            vfe = awgn_vfe(n,m,y,F,ir,jc,a,c,logz,r,sig,delta,is_array,0);
        }

        /* Update Damping */    
        if(adaptive_damp && (t > 0)){
            if(vfe > last_vfe){
                /* In this case we have observed the VFE growing */
                damp = min(damp*(1 + damp_modifier_up),damp_max);
            }else{
                /* In this case we have observed the VFE decreasing */
                damp = max(damp*(1 - damp_modifier_down),damp_min);
            }
        }

        /* Update prior parameters */
        if (learn_prior) prior((mean_removal ? n_noaux : n), r, sig, prior_prmts, a, c, NULL, 1);

        /* Update delta */
        if (learn_delta && t > 0) {
            if (!is_array) {
                delta_n = delta_d = 0; /* Sums: (w / v)^2 and (1 / v) */
                for (mu = 0; mu < (mean_removal ? m_noaux : m); mu++) {
                    delta_n += pow(w_r[mu] / v[mu], 2);
                    delta_d += (1. / v[mu]);
                }
                *delta *= (delta_n / delta_d);
            } else {
                delta_n = delta_d = 0;
                for (mu = 0; mu < (mean_removal ? m_noaux : m); mu++) {
                    delta_n += pow(w_r[mu] * delta[mu] / v[mu], 2) / delta0[mu];
                    delta_d += delta[mu] / v[mu];
                }
                gamma = delta_n / delta_d;
                if (disp) printf("\tgamma = %g\n", gamma);
                for (mu = 0; mu < m; mu++) delta[mu] = gamma * delta0[mu];
            }
        }

        /* Print some info. */
        res = 0;
        for (mu = 0; mu < m; mu++)
            res += pow(y[mu] - a_proj[mu], 2);
        res /= m;
        
        mse = 0;
        if (x) {
            if (!mean_removal) {
                for (i = 0; i < n; i++)
                    mse += pow(a[i] - x[i], 2);
                mse /= n;
            } else {
                for (i = 0; i < (n - 2); i++)
                    mse += pow(a[i] - x[i], 2);
                mse /= (n - 2);
            }
        }

        if (is_array) {
            delta_mean = 0;
            for (mu = 0; mu < m; mu++) delta_mean += delta[mu];
            delta_mean /= m;
        } else {
            delta_mean = *delta;
        }


        if(disp && adaptive_damp){
            printf("t: %3d; mse: %.4e, est noise: %.4e, rss: %.4e, diff: %.4e, damp: %.4e\n", 
                t, mse, delta_mean, res, diff / n,damp);
        }else{
            printf("t: %3d; mse: %.4e, est noise: %.4e, rss: %.4e, diff: %.4e\n", 
                t, mse, delta_mean, res, diff / n);
        }
        if (output) fprintf(output, "%d;%g;%g;%g;%g;%g;%g;%g\n", 
                t, mse, delta_mean, res, diff / n,vfe,damp,site_violation_count);
        mexEvalString("drawnow");

        if (history) {
            for (i = 0; i < n; i++) fprintf(history, "%g ", a[i]);
            fprintf(history, "\n");
        }

        /* Check for convergence */
        if (diff / n < eps) break;
    }

    /* Dealloc. structures */
    free(seq);
    free(delta0);
    free(c_proj);
    free(a_proj);
    free(g);
    free(v);
    free(w_r);
    if(calc_vfe) free(logz);
}
