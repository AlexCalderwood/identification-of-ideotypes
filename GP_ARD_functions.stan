// functions for posterior GP predictions assuming rbf kernel.
// uses analytical solution to condition on observed data
// uses ARD (automatic relevance determination, allowing different rho values
// for each parent)
functions{
  matrix rbf(vector[] x1, 
             vector[] x2, 
             real alpha, 
             vector rho) {
    
    real sq_alpha = square(alpha);
    int Nx1 = size(x1);
    int Nx2 = size(x2);
    int Nmin = min(Nx1, Nx2);
    matrix[Nx1, Nx2] K;

    for (i in 1:Nx1) {
      for (j in 1:Nx2) {
        K[i, j] = sq_alpha  
                    * exp(-0.5*dot_self((x1[i] - x2[j]) ./ rho));
      }
    }
    

    return K;
  }

  // predict y value for UNOBSERVED X_point(s)
  // all of the posterior calculations are now done within this function
	vector gp_pred_rng(vector[] x2, 
	                     vector y1, vector[] x1,
	                     real alpha, vector rho, real sigma, real delta){
    int N1 = rows(y1);
    int N2 = size(x2);
    vector[N2] f2;
    {
      matrix[N1, N1] K =  rbf(x1, x1, alpha, rho)
                         + diag_matrix(rep_vector(square(sigma), N1));
      matrix[N1, N1] L_K = cholesky_decompose(K);

      vector[N1] L_K_div_y1 = mdivide_left_tri_low(L_K, y1);
      vector[N1] K_div_y1 = mdivide_right_tri_low(L_K_div_y1', L_K)';
      matrix[N1, N2] k_x1_x2 = rbf(x1, x2, alpha, rho);
      vector[N2] f2_mu = (k_x1_x2' * K_div_y1);
      matrix[N1, N2] v_pred = mdivide_left_tri_low(L_K, k_x1_x2);
      matrix[N2, N2] cov_f2 = rbf(x2, x2, alpha, rho) - v_pred' * v_pred
                              + diag_matrix(rep_vector(delta, N2));
      f2 = multi_normal_rng(f2_mu, cov_f2);
    }
    return f2;
	}
}
