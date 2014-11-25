#include "../swamp.h"

double awgn_vfe	(size_t n, size_t m, 
			double *y,
			double *F, int *ir, int *jc,
			double *a, double *c, double *log_z,
			double *r, double *sig, 
			double *delta, int is_array){
	// This is currently really inefficient since these projections
	// are being recalculated for every coefficient. A smart person
	// would write this function such that it is a difference between
	// a current update and the new setting of variables at site "i".

	double dkl_qp = 0;
	double dkl_mout = 0;
	double *a_proj, *c_proj;
	unsigned int i, mu, idx;
	double vfe_estimate = 0;
	int sparse_mode = 1;

	if(ir==NULL){
		sparse_mode = 0;
	}

	/* Calculate what the projections of a and c are */
	a_proj = malloc(sizeof(double) * m);
    c_proj = malloc(sizeof(double) * m);
	for (mu = 0; mu < m; mu++) a_proj[mu] = c_proj[mu] = 0;
    if(sparse_mode){	    
	    for (i = 0; i < n; i++)
	        for (idx = jc[i]; idx < jc[i + 1]; idx++) {
	            a_proj[ ir[idx] ] += F[idx] * a[i];
	            c_proj[ ir[idx] ] += (F[idx] * F[idx]) * c[i];
	        }
    }else{
        for (i = 0; i < n; i++)
            for (mu = 0; mu < m; mu++) { 
                idx = mu + i * m;
                a_proj[mu] += F[idx] * a[i];
                c_proj[mu] += (F[idx] * F[idx]) * c[i];
            }    	
    }

	/* Estimate dkl_mout */
	for(mu = 0; mu < m; mu++){
		dkl_mout += (y[mu] - a_proj[mu])*(y[mu] - a_proj[mu])*(1/(2*(is_array ? delta[mu] : *delta)));
		dkl_mout += 0.5*log(2*M_PI*(is_array ? delta[mu] : *delta) + 2*M_PI*c_proj[mu]);
	}

	/* Estimate dkl_qp */
	for(i = 0; i < n; i++){
		dkl_qp += (c[i] + (a[i] - r[i])*(a[i] - r[i]))/(2*sig[i]);
		dkl_qp += log_z[i];
	}

	vfe_estimate = dkl_mout - dkl_qp;

	free(a_proj);
	free(c_proj);

	return(vfe_estimate);
}