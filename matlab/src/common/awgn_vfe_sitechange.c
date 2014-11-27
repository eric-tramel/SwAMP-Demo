#include "../swamp.h"

double awgn_vfe_sitechange	(size_t n, size_t m, 
			double *y,
			double *F, int *ir, int *jc,
			double a_site, double c_site, double log_z_site, double r_site, double sig_site, 
			double *a, double *c, double *log_z, double *r, double *sig, 
			double *delta, int is_array,
			int local_idx){
	/*------- Description -----------
	Here, we calculate the change in the VFE when a given size changes by a specified amount.
	We assume that the projections of {a} and {c} have not yet been calculated. So, the old
	values we are given should be entire vector of variables at each site from the last state.
	The new values we are given should be the *scalar* values for the site we want to update.
	--------------------------------*/

	double delta_mu_plus_c_proj;
	double delta_mu;
	double dkl_qp_change = 0;
	double dkl_mout_change = 0;
	double *a_proj, *c_proj;
	unsigned int i, mu, idx;
	double vfe_change = 0;
	int sparse_mode = 1;

	double a_diff, c_diff;
	double *a_diff_proj,*c_diff_proj;

	/* Detect Sparse Matrix Mode */
	if(ir==NULL){
		sparse_mode = 0;
	}

	/* Calculate (Fa) and (F^2c) */
	// Note: To be more efficient, one can write the AMP iteration in such a way that the values of
	//       (Fa) an (F^2c) are tracked throughout the site sweep and then could be passed to this
	//       function so as to avoid this "full" computation of the projection at each site.
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

    /* Calculate Differences */
    a_diff = a_site - a[local_idx];
    c_diff = c_site - c[local_idx];

    /* Calculate difference projections */
	a_diff_proj = malloc(sizeof(double) * m);
    c_diff_proj = malloc(sizeof(double) * m);
    if(sparse_mode){
        for (idx = jc[local_idx]; idx < jc[local_idx + 1]; idx++) {
            a_diff_proj[ ir[idx] ] = F[idx] * a_diff;
            c_diff_proj[ ir[idx] ] = (F[idx] * F[idx]) * c_diff;
        }
    }else{
	    for (mu = 0; mu < m; mu++) { 
	        idx = mu + local_idx * m;
	        a_diff_proj[mu] = F[idx] * a_diff;
	        c_diff_proj[mu] = (F[idx] * F[idx]) * c_diff;
	    }    	
    }

	/* Estimate change in dkl_mout */
	for(mu = 0; mu < m; mu++){
		delta_mu = (is_array ? delta[mu] : *delta);
		delta_mu_plus_c_proj = delta_mu + c_proj[mu];
		dkl_mout_change += (1/delta_mu)*a_diff_proj[mu]*(y[mu] - a_proj[mu]) - (1/(2*delta_mu))*a_diff_proj[mu]*a_diff_proj[mu];	
		dkl_mout_change += 0.5*log(delta_mu_plus_c_proj) - 0.5*log(delta_mu_plus_c_proj + c_diff_proj[mu]);							// Is subtracting better than dividing?
	}

	/* Estimate change in dkl_qp */
	// From the original
	dkl_qp_change += (c[local_idx] + (a[local_idx] - r[local_idx])*(a[local_idx] - r[local_idx]))/(2*sig[local_idx]);
	dkl_qp_change += log_z[local_idx];
	// Minus the new values
	dkl_qp_change -= (c_site + (a_site - r_site)*(a_site - r_site))/(2*sig_site);
	dkl_qp_change -= log_z_site;

	/* Calculate final change */
	vfe_change = dkl_mout_change - dkl_qp_change;

	free(a_proj);
	free(c_proj);
	free(a_diff_proj);
	free(c_diff_proj);

	return(vfe_change);
}