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

data{
	// Define variables in data
	//X_node = parent data of node
	//D_node = num parents of node
	int<lower=1> D_numFlowers;
	int<lower=1> D_timeToFlower;
	int<lower=1> D_ovaryLen;
	int<lower=1> D_ovuleACov;
	int<lower=1> D_ovuleNum;
	int<lower=1> D_gynLen;
	int<lower=1> D_percAbortedMain;
	int<lower=1> D_ovuleArea;
	int<lower=1> D_numPodsMain;
	int<lower=1> D_podLen;
	int<lower=1> D_percAborted2ary;
	int<lower=1> D_beakLen;
	int<lower=1> D_numPods2ary;
	int<lower=1> D_timeToMature;
	int<lower=1> D_tenPodNumber;
	int<lower=1> D_tenPodArea;
	int<lower=1> D_seedACov;
	int<lower=1> D_tenPodWeight;
	int<lower=1> D_totSeedArea;
	int<lower=1> D_totCompactness;
	int<lower=1> D_totSeedNum;
	int<lower=1> D_tenPodCompactness;
	int<lower=1> D_TGW;
	int<lower=1> D_totSeedW;
	int<lower=1> D_oilContent;
	int<lower=1> N;
	vector[N] numFlowers;
	vector[N] timeToFlower;
	vector[D_timeToFlower] X_timeToFlower[N];
	vector[N] ovaryLen;
	vector[D_ovaryLen] X_ovaryLen[N];
	vector[N] ovuleACov;
	vector[D_ovuleACov] X_ovuleACov[N];
	vector[N] ovuleNum;
	vector[N] gynLen;
	vector[N] percAbortedMain;
	vector[D_percAbortedMain] X_percAbortedMain[N];
	vector[N] ovuleArea;
	vector[D_ovuleArea] X_ovuleArea[N];
	vector[N] numPodsMain;
	vector[N] podLen;
	vector[N] percAborted2ary;
	vector[N] beakLen;
	vector[N] numPods2ary;
	vector[N] timeToMature;
	vector[D_timeToMature] X_timeToMature[N];
	vector[N] tenPodNumber;
	vector[N] tenPodArea;
	vector[N] seedACov;
	vector[N] tenPodWeight;
	vector[N] totSeedArea;
	vector[D_totSeedArea] X_totSeedArea[N];
	vector[N] totCompactness;
	vector[N] totSeedNum;
	vector[N] tenPodCompactness;
	vector[N] TGW;
	vector[N] totSeedW;
	vector[N] oilContent;
	vector[D_numFlowers] X_numFlowers[N];
	vector[D_percAborted2ary] X_percAborted2ary[N];
	vector[D_tenPodArea] X_tenPodArea[N];
	vector[D_seedACov] X_seedACov[N];
	vector[D_tenPodWeight] X_tenPodWeight[N];
	vector[D_totCompactness] X_totCompactness[N];
	vector[D_tenPodCompactness] X_tenPodCompactness[N];
	vector[D_ovuleNum] X_ovuleNum[N];
	vector[D_gynLen] X_gynLen[N];
	vector[D_beakLen] X_beakLen[N];
	vector[D_oilContent] X_oilContent[N];
	vector[D_numPodsMain] X_numPodsMain[N];
	vector[D_podLen] X_podLen[N];
	vector[D_tenPodNumber] X_tenPodNumber[N];
	vector[D_TGW] X_TGW[N];
	vector[D_numPods2ary] X_numPods2ary[N];
	vector[D_totSeedNum] X_totSeedNum[N];
	vector[D_totSeedW] X_totSeedW[N];
}

parameters {
	// Define parameters
	vector<lower=0>[D_numFlowers] Rho_numFlowers; //for automatic relevance determination, each parent gets own Rho
	real<lower=0> Alpha_numFlowers;
	real<lower=0> Sigma_numFlowers;

	vector<lower=0>[D_timeToFlower] Rho_timeToFlower; //for automatic relevance determination, each parent gets own Rho
	real<lower=0> Alpha_timeToFlower;
	real<lower=0> Sigma_timeToFlower;

	vector<lower=0>[D_ovaryLen] Rho_ovaryLen; //for automatic relevance determination, each parent gets own Rho
	real<lower=0> Alpha_ovaryLen;
	real<lower=0> Sigma_ovaryLen;

	vector<lower=0>[D_ovuleACov] Rho_ovuleACov; //for automatic relevance determination, each parent gets own Rho
	real<lower=0> Alpha_ovuleACov;
	real<lower=0> Sigma_ovuleACov;

	vector<lower=0>[D_ovuleNum] Rho_ovuleNum; //for automatic relevance determination, each parent gets own Rho
	real<lower=0> Alpha_ovuleNum;
	real<lower=0> Sigma_ovuleNum;

	vector<lower=0>[D_gynLen] Rho_gynLen; //for automatic relevance determination, each parent gets own Rho
	real<lower=0> Alpha_gynLen;
	real<lower=0> Sigma_gynLen;

	vector<lower=0>[D_percAbortedMain] Rho_percAbortedMain; //for automatic relevance determination, each parent gets own Rho
	real<lower=0> Alpha_percAbortedMain;
	real<lower=0> Sigma_percAbortedMain;

	vector<lower=0>[D_ovuleArea] Rho_ovuleArea; //for automatic relevance determination, each parent gets own Rho
	real<lower=0> Alpha_ovuleArea;
	real<lower=0> Sigma_ovuleArea;

	vector<lower=0>[D_numPodsMain] Rho_numPodsMain; //for automatic relevance determination, each parent gets own Rho
	real<lower=0> Alpha_numPodsMain;
	real<lower=0> Sigma_numPodsMain;

	vector<lower=0>[D_podLen] Rho_podLen; //for automatic relevance determination, each parent gets own Rho
	real<lower=0> Alpha_podLen;
	real<lower=0> Sigma_podLen;

	vector<lower=0>[D_percAborted2ary] Rho_percAborted2ary; //for automatic relevance determination, each parent gets own Rho
	real<lower=0> Alpha_percAborted2ary;
	real<lower=0> Sigma_percAborted2ary;

	vector<lower=0>[D_beakLen] Rho_beakLen; //for automatic relevance determination, each parent gets own Rho
	real<lower=0> Alpha_beakLen;
	real<lower=0> Sigma_beakLen;

	vector<lower=0>[D_numPods2ary] Rho_numPods2ary; //for automatic relevance determination, each parent gets own Rho
	real<lower=0> Alpha_numPods2ary;
	real<lower=0> Sigma_numPods2ary;

	vector<lower=0>[D_timeToMature] Rho_timeToMature; //for automatic relevance determination, each parent gets own Rho
	real<lower=0> Alpha_timeToMature;
	real<lower=0> Sigma_timeToMature;

	vector<lower=0>[D_tenPodNumber] Rho_tenPodNumber; //for automatic relevance determination, each parent gets own Rho
	real<lower=0> Alpha_tenPodNumber;
	real<lower=0> Sigma_tenPodNumber;

	vector<lower=0>[D_tenPodArea] Rho_tenPodArea; //for automatic relevance determination, each parent gets own Rho
	real<lower=0> Alpha_tenPodArea;
	real<lower=0> Sigma_tenPodArea;

	vector<lower=0>[D_seedACov] Rho_seedACov; //for automatic relevance determination, each parent gets own Rho
	real<lower=0> Alpha_seedACov;
	real<lower=0> Sigma_seedACov;

	vector<lower=0>[D_tenPodWeight] Rho_tenPodWeight; //for automatic relevance determination, each parent gets own Rho
	real<lower=0> Alpha_tenPodWeight;
	real<lower=0> Sigma_tenPodWeight;

	vector<lower=0>[D_totSeedArea] Rho_totSeedArea; //for automatic relevance determination, each parent gets own Rho
	real<lower=0> Alpha_totSeedArea;
	real<lower=0> Sigma_totSeedArea;

	vector<lower=0>[D_totCompactness] Rho_totCompactness; //for automatic relevance determination, each parent gets own Rho
	real<lower=0> Alpha_totCompactness;
	real<lower=0> Sigma_totCompactness;

	vector<lower=0>[D_totSeedNum] Rho_totSeedNum; //for automatic relevance determination, each parent gets own Rho
	real<lower=0> Alpha_totSeedNum;
	real<lower=0> Sigma_totSeedNum;

	vector<lower=0>[D_tenPodCompactness] Rho_tenPodCompactness; //for automatic relevance determination, each parent gets own Rho
	real<lower=0> Alpha_tenPodCompactness;
	real<lower=0> Sigma_tenPodCompactness;

	vector<lower=0>[D_TGW] Rho_TGW; //for automatic relevance determination, each parent gets own Rho
	real<lower=0> Alpha_TGW;
	real<lower=0> Sigma_TGW;

	vector<lower=0>[D_totSeedW] Rho_totSeedW; //for automatic relevance determination, each parent gets own Rho
	real<lower=0> Alpha_totSeedW;
	real<lower=0> Sigma_totSeedW;

	vector<lower=0>[D_oilContent] Rho_oilContent; //for automatic relevance determination, each parent gets own Rho
	real<lower=0> Alpha_oilContent;
	real<lower=0> Sigma_oilContent;

}

model {
	// Calculate Kernels
	matrix[N, N] K_numFlowers = rbf(X_numFlowers, X_numFlowers, Alpha_numFlowers, Rho_numFlowers)
			 + diag_matrix(rep_vector(square(Sigma_numFlowers), N));
	matrix[N, N] K_L_numFlowers = cholesky_decompose(K_numFlowers);
	matrix[N, N] K_timeToFlower = rbf(X_timeToFlower, X_timeToFlower, Alpha_timeToFlower, Rho_timeToFlower)
			 + diag_matrix(rep_vector(square(Sigma_timeToFlower), N));
	matrix[N, N] K_L_timeToFlower = cholesky_decompose(K_timeToFlower);
	matrix[N, N] K_ovaryLen = rbf(X_ovaryLen, X_ovaryLen, Alpha_ovaryLen, Rho_ovaryLen)
			 + diag_matrix(rep_vector(square(Sigma_ovaryLen), N));
	matrix[N, N] K_L_ovaryLen = cholesky_decompose(K_ovaryLen);
	matrix[N, N] K_ovuleACov = rbf(X_ovuleACov, X_ovuleACov, Alpha_ovuleACov, Rho_ovuleACov)
			 + diag_matrix(rep_vector(square(Sigma_ovuleACov), N));
	matrix[N, N] K_L_ovuleACov = cholesky_decompose(K_ovuleACov);
	matrix[N, N] K_ovuleNum = rbf(X_ovuleNum, X_ovuleNum, Alpha_ovuleNum, Rho_ovuleNum)
			 + diag_matrix(rep_vector(square(Sigma_ovuleNum), N));
	matrix[N, N] K_L_ovuleNum = cholesky_decompose(K_ovuleNum);
	matrix[N, N] K_gynLen = rbf(X_gynLen, X_gynLen, Alpha_gynLen, Rho_gynLen)
			 + diag_matrix(rep_vector(square(Sigma_gynLen), N));
	matrix[N, N] K_L_gynLen = cholesky_decompose(K_gynLen);
	matrix[N, N] K_percAbortedMain = rbf(X_percAbortedMain, X_percAbortedMain, Alpha_percAbortedMain, Rho_percAbortedMain)
			 + diag_matrix(rep_vector(square(Sigma_percAbortedMain), N));
	matrix[N, N] K_L_percAbortedMain = cholesky_decompose(K_percAbortedMain);
	matrix[N, N] K_ovuleArea = rbf(X_ovuleArea, X_ovuleArea, Alpha_ovuleArea, Rho_ovuleArea)
			 + diag_matrix(rep_vector(square(Sigma_ovuleArea), N));
	matrix[N, N] K_L_ovuleArea = cholesky_decompose(K_ovuleArea);
	matrix[N, N] K_numPodsMain = rbf(X_numPodsMain, X_numPodsMain, Alpha_numPodsMain, Rho_numPodsMain)
			 + diag_matrix(rep_vector(square(Sigma_numPodsMain), N));
	matrix[N, N] K_L_numPodsMain = cholesky_decompose(K_numPodsMain);
	matrix[N, N] K_podLen = rbf(X_podLen, X_podLen, Alpha_podLen, Rho_podLen)
			 + diag_matrix(rep_vector(square(Sigma_podLen), N));
	matrix[N, N] K_L_podLen = cholesky_decompose(K_podLen);
	matrix[N, N] K_percAborted2ary = rbf(X_percAborted2ary, X_percAborted2ary, Alpha_percAborted2ary, Rho_percAborted2ary)
			 + diag_matrix(rep_vector(square(Sigma_percAborted2ary), N));
	matrix[N, N] K_L_percAborted2ary = cholesky_decompose(K_percAborted2ary);
	matrix[N, N] K_beakLen = rbf(X_beakLen, X_beakLen, Alpha_beakLen, Rho_beakLen)
			 + diag_matrix(rep_vector(square(Sigma_beakLen), N));
	matrix[N, N] K_L_beakLen = cholesky_decompose(K_beakLen);
	matrix[N, N] K_numPods2ary = rbf(X_numPods2ary, X_numPods2ary, Alpha_numPods2ary, Rho_numPods2ary)
			 + diag_matrix(rep_vector(square(Sigma_numPods2ary), N));
	matrix[N, N] K_L_numPods2ary = cholesky_decompose(K_numPods2ary);
	matrix[N, N] K_timeToMature = rbf(X_timeToMature, X_timeToMature, Alpha_timeToMature, Rho_timeToMature)
			 + diag_matrix(rep_vector(square(Sigma_timeToMature), N));
	matrix[N, N] K_L_timeToMature = cholesky_decompose(K_timeToMature);
	matrix[N, N] K_tenPodNumber = rbf(X_tenPodNumber, X_tenPodNumber, Alpha_tenPodNumber, Rho_tenPodNumber)
			 + diag_matrix(rep_vector(square(Sigma_tenPodNumber), N));
	matrix[N, N] K_L_tenPodNumber = cholesky_decompose(K_tenPodNumber);
	matrix[N, N] K_tenPodArea = rbf(X_tenPodArea, X_tenPodArea, Alpha_tenPodArea, Rho_tenPodArea)
			 + diag_matrix(rep_vector(square(Sigma_tenPodArea), N));
	matrix[N, N] K_L_tenPodArea = cholesky_decompose(K_tenPodArea);
	matrix[N, N] K_seedACov = rbf(X_seedACov, X_seedACov, Alpha_seedACov, Rho_seedACov)
			 + diag_matrix(rep_vector(square(Sigma_seedACov), N));
	matrix[N, N] K_L_seedACov = cholesky_decompose(K_seedACov);
	matrix[N, N] K_tenPodWeight = rbf(X_tenPodWeight, X_tenPodWeight, Alpha_tenPodWeight, Rho_tenPodWeight)
			 + diag_matrix(rep_vector(square(Sigma_tenPodWeight), N));
	matrix[N, N] K_L_tenPodWeight = cholesky_decompose(K_tenPodWeight);
	matrix[N, N] K_totSeedArea = rbf(X_totSeedArea, X_totSeedArea, Alpha_totSeedArea, Rho_totSeedArea)
			 + diag_matrix(rep_vector(square(Sigma_totSeedArea), N));
	matrix[N, N] K_L_totSeedArea = cholesky_decompose(K_totSeedArea);
	matrix[N, N] K_totCompactness = rbf(X_totCompactness, X_totCompactness, Alpha_totCompactness, Rho_totCompactness)
			 + diag_matrix(rep_vector(square(Sigma_totCompactness), N));
	matrix[N, N] K_L_totCompactness = cholesky_decompose(K_totCompactness);
	matrix[N, N] K_totSeedNum = rbf(X_totSeedNum, X_totSeedNum, Alpha_totSeedNum, Rho_totSeedNum)
			 + diag_matrix(rep_vector(square(Sigma_totSeedNum), N));
	matrix[N, N] K_L_totSeedNum = cholesky_decompose(K_totSeedNum);
	matrix[N, N] K_tenPodCompactness = rbf(X_tenPodCompactness, X_tenPodCompactness, Alpha_tenPodCompactness, Rho_tenPodCompactness)
			 + diag_matrix(rep_vector(square(Sigma_tenPodCompactness), N));
	matrix[N, N] K_L_tenPodCompactness = cholesky_decompose(K_tenPodCompactness);
	matrix[N, N] K_TGW = rbf(X_TGW, X_TGW, Alpha_TGW, Rho_TGW)
			 + diag_matrix(rep_vector(square(Sigma_TGW), N));
	matrix[N, N] K_L_TGW = cholesky_decompose(K_TGW);
	matrix[N, N] K_totSeedW = rbf(X_totSeedW, X_totSeedW, Alpha_totSeedW, Rho_totSeedW)
			 + diag_matrix(rep_vector(square(Sigma_totSeedW), N));
	matrix[N, N] K_L_totSeedW = cholesky_decompose(K_totSeedW);
	matrix[N, N] K_oilContent = rbf(X_oilContent, X_oilContent, Alpha_oilContent, Rho_oilContent)
			 + diag_matrix(rep_vector(square(Sigma_oilContent), N));
	matrix[N, N] K_L_oilContent = cholesky_decompose(K_oilContent);

	// Priors act as regularisation
	Alpha_numFlowers ~ normal(0,1);
	Rho_numFlowers ~ inv_gamma(5, 5);
	Sigma_numFlowers ~ normal(0, 1);

	Alpha_timeToFlower ~ normal(0,1);
	Rho_timeToFlower ~ inv_gamma(5, 5);
	Sigma_timeToFlower ~ normal(0, 1);

	Alpha_ovaryLen ~ normal(0,1);
	Rho_ovaryLen ~ inv_gamma(5, 5);
	Sigma_ovaryLen ~ normal(0, 1);

	Alpha_ovuleACov ~ normal(0,1);
	Rho_ovuleACov ~ inv_gamma(5, 5);
	Sigma_ovuleACov ~ normal(0, 1);

	Alpha_ovuleNum ~ normal(0,1);
	Rho_ovuleNum ~ inv_gamma(5, 5);
	Sigma_ovuleNum ~ normal(0, 1);

	Alpha_gynLen ~ normal(0,1);
	Rho_gynLen ~ inv_gamma(5, 5);
	Sigma_gynLen ~ normal(0, 1);

	Alpha_percAbortedMain ~ normal(0,1);
	Rho_percAbortedMain ~ inv_gamma(5, 5);
	Sigma_percAbortedMain ~ normal(0, 1);

	Alpha_ovuleArea ~ normal(0,1);
	Rho_ovuleArea ~ inv_gamma(5, 5);
	Sigma_ovuleArea ~ normal(0, 1);

	Alpha_numPodsMain ~ normal(0,1);
	Rho_numPodsMain ~ inv_gamma(5, 5);
	Sigma_numPodsMain ~ normal(0, 1);

	Alpha_podLen ~ normal(0,1);
	Rho_podLen ~ inv_gamma(5, 5);
	Sigma_podLen ~ normal(0, 1);

	Alpha_percAborted2ary ~ normal(0,1);
	Rho_percAborted2ary ~ inv_gamma(5, 5);
	Sigma_percAborted2ary ~ normal(0, 1);

	Alpha_beakLen ~ normal(0,1);
	Rho_beakLen ~ inv_gamma(5, 5);
	Sigma_beakLen ~ normal(0, 1);

	Alpha_numPods2ary ~ normal(0,1);
	Rho_numPods2ary ~ inv_gamma(5, 5);
	Sigma_numPods2ary ~ normal(0, 1);

	Alpha_timeToMature ~ normal(0,1);
	Rho_timeToMature ~ inv_gamma(5, 5);
	Sigma_timeToMature ~ normal(0, 1);

	Alpha_tenPodNumber ~ normal(0,1);
	Rho_tenPodNumber ~ inv_gamma(5, 5);
	Sigma_tenPodNumber ~ normal(0, 1);

	Alpha_tenPodArea ~ normal(0,1);
	Rho_tenPodArea ~ inv_gamma(5, 5);
	Sigma_tenPodArea ~ normal(0, 1);

	Alpha_seedACov ~ normal(0,1);
	Rho_seedACov ~ inv_gamma(5, 5);
	Sigma_seedACov ~ normal(0, 1);

	Alpha_tenPodWeight ~ normal(0,1);
	Rho_tenPodWeight ~ inv_gamma(5, 5);
	Sigma_tenPodWeight ~ normal(0, 1);

	Alpha_totSeedArea ~ normal(0,1);
	Rho_totSeedArea ~ inv_gamma(5, 5);
	Sigma_totSeedArea ~ normal(0, 1);

	Alpha_totCompactness ~ normal(0,1);
	Rho_totCompactness ~ inv_gamma(5, 5);
	Sigma_totCompactness ~ normal(0, 1);

	Alpha_totSeedNum ~ normal(0,1);
	Rho_totSeedNum ~ inv_gamma(5, 5);
	Sigma_totSeedNum ~ normal(0, 1);

	Alpha_tenPodCompactness ~ normal(0,1);
	Rho_tenPodCompactness ~ inv_gamma(5, 5);
	Sigma_tenPodCompactness ~ normal(0, 1);

	Alpha_TGW ~ normal(0,1);
	Rho_TGW ~ inv_gamma(5, 5);
	Sigma_TGW ~ normal(0, 1);

	Alpha_totSeedW ~ normal(0,1);
	Rho_totSeedW ~ inv_gamma(5, 5);
	Sigma_totSeedW ~ normal(0, 1);

	Alpha_oilContent ~ normal(0,1);
	Rho_oilContent ~ inv_gamma(5, 5);
	Sigma_oilContent ~ normal(0, 1);


	// Likelihood section
	numFlowers ~ multi_normal_cholesky(rep_vector(0,N), K_L_numFlowers);
	timeToFlower ~ multi_normal_cholesky(rep_vector(0,N), K_L_timeToFlower);
	ovaryLen ~ multi_normal_cholesky(rep_vector(0,N), K_L_ovaryLen);
	ovuleACov ~ multi_normal_cholesky(rep_vector(0,N), K_L_ovuleACov);
	ovuleNum ~ multi_normal_cholesky(rep_vector(0,N), K_L_ovuleNum);
	gynLen ~ multi_normal_cholesky(rep_vector(0,N), K_L_gynLen);
	percAbortedMain ~ multi_normal_cholesky(rep_vector(0,N), K_L_percAbortedMain);
	ovuleArea ~ multi_normal_cholesky(rep_vector(0,N), K_L_ovuleArea);
	numPodsMain ~ multi_normal_cholesky(rep_vector(0,N), K_L_numPodsMain);
	podLen ~ multi_normal_cholesky(rep_vector(0,N), K_L_podLen);
	percAborted2ary ~ multi_normal_cholesky(rep_vector(0,N), K_L_percAborted2ary);
	beakLen ~ multi_normal_cholesky(rep_vector(0,N), K_L_beakLen);
	numPods2ary ~ multi_normal_cholesky(rep_vector(0,N), K_L_numPods2ary);
	timeToMature ~ multi_normal_cholesky(rep_vector(0,N), K_L_timeToMature);
	tenPodNumber ~ multi_normal_cholesky(rep_vector(0,N), K_L_tenPodNumber);
	tenPodArea ~ multi_normal_cholesky(rep_vector(0,N), K_L_tenPodArea);
	seedACov ~ multi_normal_cholesky(rep_vector(0,N), K_L_seedACov);
	tenPodWeight ~ multi_normal_cholesky(rep_vector(0,N), K_L_tenPodWeight);
	totSeedArea ~ multi_normal_cholesky(rep_vector(0,N), K_L_totSeedArea);
	totCompactness ~ multi_normal_cholesky(rep_vector(0,N), K_L_totCompactness);
	totSeedNum ~ multi_normal_cholesky(rep_vector(0,N), K_L_totSeedNum);
	tenPodCompactness ~ multi_normal_cholesky(rep_vector(0,N), K_L_tenPodCompactness);
	TGW ~ multi_normal_cholesky(rep_vector(0,N), K_L_TGW);
	totSeedW ~ multi_normal_cholesky(rep_vector(0,N), K_L_totSeedW);
	oilContent ~ multi_normal_cholesky(rep_vector(0,N), K_L_oilContent);
}

