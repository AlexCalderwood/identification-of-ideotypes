functions{
  // can't use cov_exp_quad, because want to use vector of rho. BUT maybe could 
  // get away with one only due to alpha?
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
    
    // // try and fix numerical stability issue
    // for (i in 1:Nmin) {
    //   K[i,i] = sq_alpha + delta;
    // }
    
    return K;
  }
  
	// get the mean prediction at point X2 conditioned on y1, x1
	// (not a sample from the point)
	vector gp_pred_mean(vector[] x2, 
	                     vector y1, vector[] x1,
	                     real alpha, vector rho, real sigma, real delta){
    int N1 = rows(y1);
    int N2 = size(x2);
    vector[N2] f2_mu;
    {
      //matrix[N1, N1] K =  cov_exp_quad(x1, alpha, rho)
      matrix[N1, N1] K =  rbf(x1, x1, alpha, rho)
                         + diag_matrix(rep_vector(square(sigma), N1));
      matrix[N1, N1] L_K = cholesky_decompose(K);

      vector[N1] L_K_div_y1 = mdivide_left_tri_low(L_K, y1);
      vector[N1] K_div_y1 = mdivide_right_tri_low(L_K_div_y1', L_K)';
      //matrix[N1, N2] k_x1_x2 = cov_exp_quad(x1, x2, alpha, rho);
      matrix[N1, N2] k_x1_x2 = rbf(x1, x2, alpha, rho);
      // mean of predsd
      f2_mu = (k_x1_x2' * K_div_y1);
    }
    return f2_mu;
	}
	
	 // calculate cov matrix of preds at locations x2
		matrix gp_pred_cov(vector[] x2, 
	                     vector y1, vector[] x1,
	                     real alpha, vector rho, real sigma, real delta){
    int N1 = rows(y1);
    int N2 = size(x2);
    matrix[N2, N2] cov_f2;
    {
      //matrix[N1, N1] K =  cov_exp_quad(x1, alpha, rho)
      matrix[N1, N1] K =  rbf(x1, x1, alpha, rho)
                         + diag_matrix(rep_vector(square(sigma), N1));
      matrix[N1, N1] L_K = cholesky_decompose(K);

      vector[N1] L_K_div_y1 = mdivide_left_tri_low(L_K, y1);
      vector[N1] K_div_y1 = mdivide_right_tri_low(L_K_div_y1', L_K)';
      //matrix[N1, N2] k_x1_x2 = cov_exp_quad(x1, x2, alpha, rho);
      matrix[N1, N2] k_x1_x2 = rbf(x1, x2, alpha, rho);
      vector[N2] f2_mu = (k_x1_x2' * K_div_y1);
      matrix[N1, N2] v_pred = mdivide_left_tri_low(L_K, k_x1_x2);
      //matrix[N2, N2] cov_f2 =   cov_exp_quad(x2, x2, alpha, rho) - v_pred' * v_pred
      
      cov_f2 = rbf(x2, x2, alpha, rho) - v_pred' * v_pred
                              + diag_matrix(rep_vector(delta, N2));
    }
    return cov_f2;
	}
}

data {
  int<lower=1> N; // num observations
  int<lower=1> D_totSeedW; // number of covariates
  vector[D_totSeedW] X_totSeedW[N];
  vector[N] totSeedW;
  
  int<lower=1> N_pred; // num predictions to make
  vector[D_totSeedW] X_pred[N_pred];
  
  vector<lower=0>[D_totSeedW] Rho_totSeedW;
  real<lower=0> Alpha_totSeedW;
  real<lower=0> Sigma_totSeedW;
}

parameters{}

model{}

generated quantities {
  // predictions at points with observations
  
  vector[N] y = gp_pred_mean(X_totSeedW, totSeedW, X_totSeedW,
    Alpha_totSeedW, Rho_totSeedW, Sigma_totSeedW, 1e-10);

  // prediction at new point
  vector[N_pred] f_pred = gp_pred_mean(X_pred, totSeedW, X_totSeedW, 
  Alpha_totSeedW, Rho_totSeedW, Sigma_totSeedW, 1e-10);
  
  matrix[N_pred, N_pred] cov_pred = gp_pred_cov(X_pred, totSeedW, X_totSeedW, 
  Alpha_totSeedW, Rho_totSeedW, Sigma_totSeedW, 1e-10);

  
  // calculate EI, following 
  // http://krasserm.github.io/2018/03/21/bayesian-optimization/
  // vector[N_pred] ei;
  // {
  // real y_opt = max(y);
  // vector[N_pred] imp = f_pred - y_opt - epsilon;
  // vector[N_pred] Z = imp / Sigma_totSeedW;
  // 
  // vector[N_pred] zeros = rep_vector(0, N_pred);
  // vector[N_pred] ones = rep_vector(1, N_pred);
  // ei = normal_cdf(Z, zeros, ones); 
  // (Sigma_totSeedW * exp(normal_lpdf(Z | 0, 1)));
  // }
}







