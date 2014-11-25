#include "cswgamp.h"

/* Interface b/w MATLAB and C */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
    double *ry, *rF, *cy, *cF;          /* Input/output related */
    double *ch_prmts, *rpr_prmts, *cpr_prmts, *rinit, *cinit, t_max, eps, damp; 
    double *ra, *ca, *c, *rr, *cr, *sig, *diff;
    complex double *y, *F, *a, *r;
    complex double *pr_prmts;
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
    ry = mxGetPr(prhs[0]); cy = mxGetPi(prhs[0]);
    rF = mxGetPr(prhs[1]); cF = mxGetPi(prhs[1]);
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
            opt_ou = opt_hi = opt_si = NULL;
    }

    if (opt_ch) {                       /* Channel type */
        if (strcmp(mxArrayToString(opt_ch), "cpr") == 0) 
            channel = channel_cpr;
        else
            channel = channel_cgaussian;
    }

    if (opt_cp) {                       /* Channel parameters */
        plhs[4] = mxDuplicateArray(opt_cp);
        ch_prmts = mxGetPr(plhs[4]);
        learn_channel = opt_lc ? *mxGetPr(opt_lc) : 0;
    } else {
        channel = channel_cgaussian;
        plhs[4] = mxCreateDoubleMatrix(2, 1, mxREAL);
        ch_prmts = mxGetPr(plhs[4]);
        ch_prmts[0] = 1.;
        learn_channel = 1;
    }

    prior = prior_cgb;                 /* Prior */
    pr_prmts = malloc(sizeof(complex double) * 3);
    if (opt_pm) {
        plhs[5] = mxDuplicateArray(opt_pm);
        rpr_prmts = mxGetPr(plhs[5]); cpr_prmts = mxGetPi(plhs[5]);
        for (key = 0; key < 3; key++)
            pr_prmts[key] = rpr_prmts[key] + (cpr_prmts ? cpr_prmts[key] : 0);
        learn_prior = opt_lp ? *mxGetPr(opt_lp) : 0;
    } else {
        prior = prior_cgb;
        plhs[5] = mxCreateDoubleMatrix(3, 1, mxCOMPLEX);
        rpr_prmts = mxGetPr(plhs[5]); cpr_prmts = mxGetPi(plhs[5]);
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
    plhs[0] = mxCreateDoubleMatrix(n, 1, mxCOMPLEX);
    plhs[1] = mxCreateDoubleMatrix(n, 1, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(n, 1, mxCOMPLEX);
    plhs[3] = mxCreateDoubleMatrix(n, 1, mxREAL);
    plhs[6] = mxCreateDoubleMatrix(1, 1, mxREAL);

    ra = mxGetPr(plhs[0]); ca = mxGetPi(plhs[0]);
    c = mxGetPr(plhs[1]);
    rr = mxGetPr(plhs[2]); cr = mxGetPi(plhs[2]);
    sig = mxGetPr(plhs[3]);
    diff = mxGetPr(plhs[6]);

    if (opt_in) {
        rinit = mxGetPr(opt_in); cinit = mxGetPi(opt_in);
        for (i = 0; i < n; i++) {
            ra[i] = rinit[i], c[i] = rinit[n + i];
            ca[i] = cinit ? cinit[i] : 0.;
        }
    } else {
        for (i = 0; i < n; i++) ra[i] = 0., ca[i] = 0., c[i] = 1.;
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

        /* Alloc. and init. complex structures */
        y = mxMalloc(sizeof(complex double) * m);
        F = mxMalloc(sizeof(complex double) * nnz);
        a = mxMalloc(sizeof(complex double) * n);
        r = mxMalloc(sizeof(complex double) * n);

        for (mu = 0; mu < m; mu++) y[mu] = ry[mu] + (cy ? cy[mu] * I : 0);
        for (key = 0; key < nnz; key++) F[key] = rF[key] + (cF ? cF[key] * I : 0);
        for (i = 0; i < n; i++) a[i] = ra[i] + (ca ? ca[i] * I : 0);
        for (i = 0; i < n; i++) r[i] = rr[i] + (cr ? cr[i] * I : 0);

        cgamp(n, m, y, F, ir, jc, 
            channel, ch_prmts, learn_channel,
            prior, pr_prmts, learn_prior, 
            t_max, eps, damp, disp, output, history, x,
            a, c, r, sig, diff);

        for (i = 0; i < n; i++) ra[i] = creal(a[i]), rr[i] = creal(r[i]);
        if (ca && cr)
            for (i = 0; i < n; i++) ca[i] = cimag(a[i]), cr[i] = cimag(r[i]);

        for (key = 0; key < 3; key++) rpr_prmts[key] = creal(pr_prmts[key]);
        if (cpr_prmts)
            for (key = 0; key < 3; key++) cpr_prmts[key] = cimag(pr_prmts[key]);
    } else {
        mexErrMsgTxt("Convert F to sparse before using this solver.");
    };

    /* Dealloc. structures */
    mxFree(r);
    mxFree(a);
    mxFree(F);
    mxFree(y);
    mxFree(jc);
    mxFree(ir);

    if (output) fclose(output);
}
