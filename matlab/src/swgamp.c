#include "swgamp.h"

/* Interface b/w MATLAB and C */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
    double *y, *F;                      /* Input/output related */
    double *ch_prmts, *pr_prmts, *init, t_max, eps, damp; 
    double *a, *c, *r, *sig;
    int learn_channel, learn_prior;

    void (*channel), (*prior);          /* Auxiliary variables */
    int *ir, *jc;

    mwIndex *ir_, *jc_;
    mxArray *opt_tm, *opt_ep, *opt_ch, *opt_cp, *opt_lc,
            *opt_pr, *opt_pm, *opt_lp, *opt_in, *opt_da, 
            *opt_di, *opt_ou, *opt_hi, *opt_si;

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
        opt_ch = mxGetField(prhs[2], 0, "channelType");
        opt_cp = mxGetField(prhs[2], 0, "channelPrmts");
        opt_lc = mxGetField(prhs[2], 0, "learnChannel");
        opt_pr = mxGetField(prhs[2], 0, "priorDistr");
        opt_pm = mxGetField(prhs[2], 0, "priorPrmts");
        opt_lp = mxGetField(prhs[2], 0, "learnPrior");
        opt_tm = mxGetField(prhs[2], 0, "maxIter");
        opt_ep = mxGetField(prhs[2], 0, "prec");
        opt_in = mxGetField(prhs[2], 0, "initState");
        opt_da = mxGetField(prhs[2], 0, "damp");
        opt_di = mxGetField(prhs[2], 0, "display");
        opt_ou = mxGetField(prhs[2], 0, "output");
        opt_hi = mxGetField(prhs[2], 0, "history");
        opt_si = mxGetField(prhs[2], 0, "signal");
    } else {
        opt_tm = opt_ep = opt_ch = opt_cp = opt_lc =
            opt_pr = opt_pm = opt_lp = opt_in = opt_da = 
            opt_di = opt_ou = opt_hi = NULL;
    }

    if (opt_ch) {                       /* Channel type */
        if (strcmp(mxArrayToString(opt_ch), "probit") == 0) 
            channel = channel_probit;
        else
            channel = channel_gaussian;
    }

    if (opt_cp) {                       /* Channel parameters */
        plhs[4] = mxDuplicateArray(opt_cp);
        ch_prmts = mxGetPr(plhs[4]);
        learn_channel = opt_lc ? *mxGetPr(opt_lc) : 0;
    } else if (strcmp(mxArrayToString(opt_ch), "bit") != 0) {
            channel = channel_gaussian;
            plhs[4] = mxCreateDoubleMatrix(2, 1, mxREAL);
            ch_prmts = mxGetPr(plhs[4]);
            ch_prmts[0] = 1.;
            learn_channel = 1;
    }

    if (opt_pr) {                       /* Prior distribution */
        if (strcmp(mxArrayToString(opt_pr), "binary") == 0)
            prior = prior_binary;
        else
            prior = prior_gb;
    } else {
        prior = prior_gb; 
    }

    if (opt_pm) {                       /* Prior parameters */
        plhs[5] = mxDuplicateArray(opt_pm);
        pr_prmts = mxGetPr(plhs[5]);
        learn_prior = opt_lp ? *mxGetPr(opt_lp) : 0;
    } else {
        prior = prior_gb;
        plhs[5] = mxCreateDoubleMatrix(3, 1, mxREAL);
        pr_prmts = mxGetPr(plhs[5]);
        pr_prmts[0] = 0.5, pr_prmts[1] = 0.0, pr_prmts[2] = 1.0;
        learn_prior = 1; 
    }

    t_max = opt_tm ? *mxGetPr(opt_tm) : 250;
    eps = opt_ep ? *mxGetPr(opt_ep) : 1e-13;
    damp = opt_da ? *mxGetPr(opt_da) : 0;
    disp = opt_di ? *mxGetPr(opt_di) : 1;
    output = opt_ou ? fopen(mxArrayToString(opt_ou), "w") : NULL;
    history = opt_hi ? fopen(mxArrayToString(opt_hi), "w") : NULL;
    x = opt_si ? mxGetPr(mxDuplicateArray(opt_si)) : NULL;

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

        gamp(n, m, y, F, ir, jc, 
            channel, ch_prmts, learn_channel,
            prior, pr_prmts, learn_prior, 
            t_max, eps, damp, disp, output, history, x,
            a, c, r, sig);
    } else {
        mexErrMsgTxt("Convert F to sparse before using this solver.");
    };

    /* Dealloc. structures */
    if(output) fclose(output);
}



