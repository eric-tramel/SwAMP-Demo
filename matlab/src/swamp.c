#include "swamp.h"

/* Interface b/w MATLAB and C */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
    double *y, *F;                      /* Input/output related */
    double *delta, *prmts;
    double *init, t_max, eps, damp; 
    double *a, *c, *r, *sig;
    int learn_delta, learn_prior;

    /* Extra Feature Flags*/
    int mean_removal;                   /* Use mean removal via auxiliary variables */
    int adaptive_damp;                  /* Use VFE to change damping per sweep */
    int calc_vfe;                       /* Save/Output per sweep VFE calc */
    int no_violations;                  /* Force VFE to decrease, strictly. */
    int site_rejection;                 /* Check for VFE violations at the site level */

    void (*prior);                      /* Auxiliary variables */
    int *ir, *jc;
    int is_array;

    mwIndex *ir_, *jc_;
    mxArray *opt_so, *opt_tm, *opt_ep, *opt_pr, *opt_pp, 
            *opt_lp, *opt_cp, *opt_lc, *opt_in, *opt_da, 
            *opt_di, *opt_ou, *opt_hi, *opt_si, *opt_mr,
            *opt_ad, *opt_cv, *opt_nv, *opt_sr;

    size_t m, n, nnz;
    unsigned int mu, i, key;

    FILE *output, *history;
    double *x;
    int disp;

    /* Handle input */
    /* TODO: check number of arguments (types and lengths as well?) */
    y = mxGetPr(prhs[0]);
    F = mxGetPr(prhs[1]);
    m = mxGetM(prhs[1]);
    n = mxGetN(prhs[1]);

    if (nrhs > 2) {
        opt_so = mxGetField(prhs[2], 0, "solver");
        opt_cp = mxGetField(prhs[2], 0, "delta");
        opt_lc = mxGetField(prhs[2], 0, "learnDelta");
        opt_pr = mxGetField(prhs[2], 0, "priorDistr");
        opt_pp = mxGetField(prhs[2], 0, "priorPrmts");
        opt_lp = mxGetField(prhs[2], 0, "learnPrior");
        opt_tm = mxGetField(prhs[2], 0, "maxIter");
        opt_ep = mxGetField(prhs[2], 0, "prec");
        opt_in = mxGetField(prhs[2], 0, "initState");
        opt_da = mxGetField(prhs[2], 0, "damp");
        opt_di = mxGetField(prhs[2], 0, "display");
        opt_ou = mxGetField(prhs[2], 0, "output");
        opt_hi = mxGetField(prhs[2], 0, "history");
        opt_si = mxGetField(prhs[2], 0, "signal");
        opt_mr = mxGetField(prhs[2], 0, "mean_removal");
        opt_ad = mxGetField(prhs[2], 0, "adaptive_damp");
        opt_cv = mxGetField(prhs[2], 0, "calc_vfe");
        opt_nv = mxGetField(prhs[2], 0, "no_violations");
        opt_sr = mxGetField(prhs[2], 0, "site_rejection");
    } else {
            opt_so = opt_tm = opt_ep = opt_pr = opt_pp = 
            opt_lp = opt_cp = opt_lc = opt_in = opt_da = 
            opt_di = opt_ou = opt_hi = opt_si = opt_mr =
            opt_ad = opt_cv = opt_nv = opt_sr = NULL;
    }

    if (opt_cp) {                       /* Delta */
        delta = mxGetPr(mxDuplicateArray(opt_cp));
        is_array = mxGetM(opt_cp) > 1;
        learn_delta = opt_lc ? *mxGetPr(opt_lc) : 0;
        if (learn_delta && is_array) learn_delta = 2;
    } else {
        delta = mxMalloc(sizeof(double));
        *delta = 1.0;
        is_array = 0;
        learn_delta = 1;
    }

    if (opt_pr) {                       /* Prior distribution */
        if (strcmp(mxArrayToString(opt_pr), "binary") == 0)
            prior = prior_binary;
        else
            prior = prior_gb;
    } else {
        prior = prior_gb; 
    }

    if (opt_pp) {                       /* Prior dist. parameters */
        plhs[4] = mxDuplicateArray(opt_pp);
        prmts = mxGetPr(plhs[4]);
        learn_prior = opt_lp ? *mxGetPr(opt_lp) : 0;
    } else {
        prior = prior_gb;
        plhs[4] = mxCreateDoubleMatrix(3, 1, mxREAL);
        prmts = mxGetPr(plhs[4]);
        prmts[0] = 0.5, prmts[1] = 0.0, prmts[2] = 1.0;
        learn_prior = 1; 
    }

    t_max = opt_tm ? *mxGetPr(opt_tm) : 250;
    eps = opt_ep ? *mxGetPr(opt_ep) : 1e-13;
    damp = opt_da ? *mxGetPr(opt_da) : 0;
    disp = opt_di ? *mxGetPr(opt_di) : 1;
    output = opt_ou ? fopen(mxArrayToString(opt_ou), "w") : NULL;
    history = opt_hi ? fopen(mxArrayToString(opt_hi), "w") : NULL;
    x = opt_si ? mxGetPr(mxDuplicateArray(opt_si)) : NULL;

    /* Extra Features */
    mean_removal = opt_mr ? *mxGetPr(opt_mr) : 0;
    adaptive_damp = opt_ad ? *mxGetPr(opt_ad) : 0;
    calc_vfe = opt_cv ? *mxGetPr(opt_cv) : 0;
    no_violations = opt_nv ? *mxGetPr(opt_nv) : 0;
    site_rejection = opt_sr ? *mxGetPr(opt_sr) : 0;
        /* Force calc_vfe flag if we need it. */
        if(adaptive_damp && ~calc_vfe) calc_vfe = 1;
        if(site_rejection && ~calc_vfe) calc_vfe = 1;        

    if(disp){
        printf("\n");
        printf("[ * Mean Removal   : %d]\n",mean_removal);
        printf("[ * Adaptive Damp  : %d]\n",adaptive_damp);
        printf("[ * Calc VFE       : %d]\n",calc_vfe);
        printf("[ * No Violations  : %d]\n",no_violations);
        printf("[ * Site Rejection : %d]\n",site_rejection);
    }


    /* Set-up output */
    plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(n, 1, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(n, 1, mxREAL);
    plhs[3] = mxCreateDoubleMatrix(n, 1, mxREAL);

    a = mxGetPr(plhs[0]);
    c = mxGetPr(plhs[1]);
    r = mxGetPr(plhs[2]);
    sig = mxGetPr(plhs[3]);

    if (opt_in) {
        init = mxGetPr(opt_in);
        for (i = 0; i < n; i++) a[i] = init[i], c[i] = init[n + i];
    } else {
        for (i = 0; i < n; i++) a[i] = 0., c[i] = 1.;
    }

    /* Run algorithm */
    if (mxIsSparse(prhs[1])) {
        /* Generate arrays w/ indexes */
        ir_ = mxGetIr(prhs[1]);
        jc_ = mxGetJc(prhs[1]);
        nnz = jc_[n];

        /* Cast indexes from mwArray to int */
        ir = mxMalloc(sizeof(int) * nnz);
        jc = mxMalloc(sizeof(int) * (n + 1));
        for (key = 0; key < nnz; key++) ir[key] = ir_[key];
        for (i = 0; i < n + 1; i++) jc[i] = jc_[i];

        if (opt_so)
            if (strcmp(mxArrayToString(opt_so), "amp_alt") == 0)
                amp_alt(n, m, y, F, ir, jc, 
                    delta, is_array, learn_delta, 
                    prior, prmts, learn_prior, 
                    t_max, eps, damp, disp, output, history, x,
                    a, c, r, sig,
                    mean_removal,calc_vfe,adaptive_damp,no_violations,site_rejection);
            else
                amp(n, m, y, F, ir, jc, 
                    delta, is_array, learn_delta, 
                    prior, prmts, learn_prior, 
                    t_max, eps, damp, disp, output, history, x,
                    a, c, r, sig,
                    mean_removal,calc_vfe,adaptive_damp,no_violations,site_rejection);
        else
            amp(n, m, y, F, ir, jc, 
                delta, is_array, learn_delta, 
                prior, prmts, learn_prior, 
                t_max, eps, damp, disp, output, history, x,
                a, c, r, sig,
                mean_removal,calc_vfe,adaptive_damp,no_violations,site_rejection);
    } else {
        if (opt_so)
            if (strcmp(mxArrayToString(opt_so), "nmf") == 0)
                mexErrMsgTxt("Convert F to sparse before using this solver.");
            else
                amp_dense(n, m, y, F,
                    delta, is_array, learn_delta, 
                    prior, prmts, learn_prior, 
                    t_max, eps, damp, disp, output, history, x,
                    a, c, r, sig,
                    mean_removal,calc_vfe,adaptive_damp,no_violations,site_rejection);
        else
            amp_dense(n, m, y, F,
                delta, is_array, learn_delta, 
                prior, prmts, learn_prior, 
                t_max, eps, damp, disp, output, history, x,
                a, c, r, sig,
                mean_removal,calc_vfe,adaptive_damp,no_violations,site_rejection);
    }

    /* Dealloc. structures */
    if (output) fclose(output);
}
