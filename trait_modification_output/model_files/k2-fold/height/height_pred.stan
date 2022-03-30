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
	int<lower=1> NPred;
	real<lower=0> Alpha_numFlowers;
	real<lower=0> Sigma_numFlowers;
	vector<lower=0>[D_timeToFlower] Rho_timeToFlower;
	real<lower=0> Alpha_timeToFlower;
	real<lower=0> Sigma_timeToFlower;
	vector<lower=0>[D_ovaryLen] Rho_ovaryLen;
	real<lower=0> Alpha_ovaryLen;
	real<lower=0> Sigma_ovaryLen;
	vector<lower=0>[D_ovuleACov] Rho_ovuleACov;
	real<lower=0> Alpha_ovuleACov;
	real<lower=0> Sigma_ovuleACov;
	real<lower=0> Alpha_ovuleNum;
	real<lower=0> Sigma_ovuleNum;
	real<lower=0> Alpha_gynLen;
	real<lower=0> Sigma_gynLen;
	vector<lower=0>[D_percAbortedMain] Rho_percAbortedMain;
	real<lower=0> Alpha_percAbortedMain;
	real<lower=0> Sigma_percAbortedMain;
	vector<lower=0>[D_ovuleArea] Rho_ovuleArea;
	real<lower=0> Alpha_ovuleArea;
	real<lower=0> Sigma_ovuleArea;
	real<lower=0> Alpha_numPodsMain;
	real<lower=0> Sigma_numPodsMain;
	real<lower=0> Alpha_podLen;
	real<lower=0> Sigma_podLen;
	real<lower=0> Alpha_percAborted2ary;
	real<lower=0> Sigma_percAborted2ary;
	real<lower=0> Alpha_beakLen;
	real<lower=0> Sigma_beakLen;
	real<lower=0> Alpha_numPods2ary;
	real<lower=0> Sigma_numPods2ary;
	vector<lower=0>[D_timeToMature] Rho_timeToMature;
	real<lower=0> Alpha_timeToMature;
	real<lower=0> Sigma_timeToMature;
	real<lower=0> Alpha_tenPodNumber;
	real<lower=0> Sigma_tenPodNumber;
	real<lower=0> Alpha_tenPodArea;
	real<lower=0> Sigma_tenPodArea;
	real<lower=0> Alpha_seedACov;
	real<lower=0> Sigma_seedACov;
	real<lower=0> Alpha_tenPodWeight;
	real<lower=0> Sigma_tenPodWeight;
	vector<lower=0>[D_totSeedArea] Rho_totSeedArea;
	real<lower=0> Alpha_totSeedArea;
	real<lower=0> Sigma_totSeedArea;
	real<lower=0> Alpha_totCompactness;
	real<lower=0> Sigma_totCompactness;
	real<lower=0> Alpha_totSeedNum;
	real<lower=0> Sigma_totSeedNum;
	real<lower=0> Alpha_tenPodCompactness;
	real<lower=0> Sigma_tenPodCompactness;
	real<lower=0> Alpha_TGW;
	real<lower=0> Sigma_TGW;
	real<lower=0> Alpha_totSeedW;
	real<lower=0> Sigma_totSeedW;
	real<lower=0> Alpha_oilContent;
	real<lower=0> Sigma_oilContent;
	vector<lower=0>[D_numFlowers] Rho_numFlowers;
	vector<lower=0>[D_percAborted2ary] Rho_percAborted2ary;
	vector<lower=0>[D_tenPodArea] Rho_tenPodArea;
	vector<lower=0>[D_seedACov] Rho_seedACov;
	vector<lower=0>[D_tenPodWeight] Rho_tenPodWeight;
	vector<lower=0>[D_totCompactness] Rho_totCompactness;
	vector<lower=0>[D_tenPodCompactness] Rho_tenPodCompactness;
	vector<lower=0>[D_ovuleNum] Rho_ovuleNum;
	vector<lower=0>[D_gynLen] Rho_gynLen;
	vector<lower=0>[D_beakLen] Rho_beakLen;
	vector<lower=0>[D_oilContent] Rho_oilContent;
	vector<lower=0>[D_numPodsMain] Rho_numPodsMain;
	vector<lower=0>[D_podLen] Rho_podLen;
	vector<lower=0>[D_tenPodNumber] Rho_tenPodNumber;
	vector<lower=0>[D_TGW] Rho_TGW;
	vector<lower=0>[D_numPods2ary] Rho_numPods2ary;
	vector<lower=0>[D_totSeedNum] Rho_totSeedNum;
	vector<lower=0>[D_totSeedW] Rho_totSeedW;
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
	vector[NPred] num2aryBranch_4pred;
	vector[NPred] styleLen_4pred;
	vector[NPred] ovuleACov_4pred;
	vector[NPred] height_4pred;
}

transformed data {
	real delta = 1e-9; // needed to ensure K positive semidef
}

parameters{
}

model{
}

generated quantities {
	vector[N] numFlowers_rep;
	vector[N] timeToFlower_rep;
	vector[N] ovaryLen_rep;
	vector[N] ovuleACov_rep;
	vector[N] ovuleNum_rep;
	vector[N] gynLen_rep;
	vector[N] percAbortedMain_rep;
	vector[N] ovuleArea_rep;
	vector[N] numPodsMain_rep;
	vector[N] podLen_rep;
	vector[N] percAborted2ary_rep;
	vector[N] beakLen_rep;
	vector[N] numPods2ary_rep;
	vector[N] timeToMature_rep;
	vector[N] tenPodNumber_rep;
	vector[N] tenPodArea_rep;
	vector[N] seedACov_rep;
	vector[N] tenPodWeight_rep;
	vector[N] totSeedArea_rep;
	vector[N] totCompactness_rep;
	vector[N] totSeedNum_rep;
	vector[N] tenPodCompactness_rep;
	vector[N] TGW_rep;
	vector[N] totSeedW_rep;
	vector[N] oilContent_rep;

	vector[NPred] numFlowers_pred;
	vector[NPred] timeToFlower_pred;
	vector[NPred] ovaryLen_pred;
	vector[NPred] ovuleNum_pred;
	vector[NPred] gynLen_pred;
	vector[NPred] percAbortedMain_pred;
	vector[NPred] ovuleArea_pred;
	vector[NPred] numPodsMain_pred;
	vector[NPred] podLen_pred;
	vector[NPred] percAborted2ary_pred;
	vector[NPred] beakLen_pred;
	vector[NPred] numPods2ary_pred;
	vector[NPred] timeToMature_pred;
	vector[NPred] tenPodNumber_pred;
	vector[NPred] tenPodArea_pred;
	vector[NPred] seedACov_pred;
	vector[NPred] tenPodWeight_pred;
	vector[NPred] totSeedArea_pred;
	vector[NPred] totCompactness_pred;
	vector[NPred] totSeedNum_pred;
	vector[NPred] tenPodCompactness_pred;
	vector[NPred] TGW_pred;
	vector[NPred] totSeedW_pred;
	vector[NPred] oilContent_pred;


	// Predictions for the observed data
	{
		vector[N] f_ovuleACov_rep;

		f_ovuleACov_rep = gp_pred_rng(X_ovuleACov, ovuleACov, X_ovuleACov,
					Alpha_ovuleACov, Rho_ovuleACov, Sigma_ovuleACov, delta);

		for (n in 1:N) {
			ovuleACov_rep[n] = normal_rng(f_ovuleACov_rep[n], Sigma_ovuleACov);
		}
	}


	//Predictions for the perturbed data
	{
		// numFlowers
		vector[NPred+N] f_numFlowers_pred;
		vector[NPred+N] p_numFlowers_pred;
		vector[D_numFlowers] X_numFlowersPred[NPred]; // declare X array for prediction
		vector[D_numFlowers] X_numFlowersAll[NPred+N]; // declare X array for prediction (_pred and _rep)
		// fill in X array
		for (n in 1:NPred) {
			X_numFlowersPred[n, 1] = height_4pred[n]; //parents ordered alphabetically
			X_numFlowersPred[n, 2] = num2aryBranch_4pred[n]; //parents ordered alphabetically
		}

		X_numFlowersAll = append_array(X_numFlowers, X_numFlowersPred);

		f_numFlowers_pred = gp_pred_rng(X_numFlowersAll,numFlowers, X_numFlowers,
					Alpha_numFlowers, Rho_numFlowers, Sigma_numFlowers, delta);

		for (n in 1:N+NPred) {
			p_numFlowers_pred[n] = normal_rng(f_numFlowers_pred[n], Sigma_numFlowers);
		}

		//Split the _rep and _pred predictions
		numFlowers_rep = p_numFlowers_pred[1:N];
		numFlowers_pred = p_numFlowers_pred[N+1:N+NPred];
	}

	{
		// timeToFlower
		vector[NPred+N] f_timeToFlower_pred;
		vector[NPred+N] p_timeToFlower_pred;
		vector[D_timeToFlower] X_timeToFlowerPred[NPred]; // declare X array for prediction
		vector[D_timeToFlower] X_timeToFlowerAll[NPred+N]; // declare X array for prediction (_pred and _rep)
		// fill in X array
		for (n in 1:NPred) {
			X_timeToFlowerPred[n, 1] = height_4pred[n]; //parents ordered alphabetically
		}

		X_timeToFlowerAll = append_array(X_timeToFlower, X_timeToFlowerPred);

		f_timeToFlower_pred = gp_pred_rng(X_timeToFlowerAll,timeToFlower, X_timeToFlower,
					Alpha_timeToFlower, Rho_timeToFlower, Sigma_timeToFlower, delta);

		for (n in 1:N+NPred) {
			p_timeToFlower_pred[n] = normal_rng(f_timeToFlower_pred[n], Sigma_timeToFlower);
		}

		//Split the _rep and _pred predictions
		timeToFlower_rep = p_timeToFlower_pred[1:N];
		timeToFlower_pred = p_timeToFlower_pred[N+1:N+NPred];
	}

	{
		// ovaryLen
		vector[NPred+N] f_ovaryLen_pred;
		vector[NPred+N] p_ovaryLen_pred;
		vector[D_ovaryLen] X_ovaryLenPred[NPred]; // declare X array for prediction
		vector[D_ovaryLen] X_ovaryLenAll[NPred+N]; // declare X array for prediction (_pred and _rep)
		// fill in X array
		for (n in 1:NPred) {
			X_ovaryLenPred[n, 1] = height_4pred[n]; //parents ordered alphabetically
		}

		X_ovaryLenAll = append_array(X_ovaryLen, X_ovaryLenPred);

		f_ovaryLen_pred = gp_pred_rng(X_ovaryLenAll,ovaryLen, X_ovaryLen,
					Alpha_ovaryLen, Rho_ovaryLen, Sigma_ovaryLen, delta);

		for (n in 1:N+NPred) {
			p_ovaryLen_pred[n] = normal_rng(f_ovaryLen_pred[n], Sigma_ovaryLen);
		}

		//Split the _rep and _pred predictions
		ovaryLen_rep = p_ovaryLen_pred[1:N];
		ovaryLen_pred = p_ovaryLen_pred[N+1:N+NPred];
	}

	{
		// ovuleNum
		vector[NPred+N] f_ovuleNum_pred;
		vector[NPred+N] p_ovuleNum_pred;
		vector[D_ovuleNum] X_ovuleNumPred[NPred]; // declare X array for prediction
		vector[D_ovuleNum] X_ovuleNumAll[NPred+N]; // declare X array for prediction (_pred and _rep)
		// fill in X array
		for (n in 1:NPred) {
			X_ovuleNumPred[n, 1] = num2aryBranch_4pred[n]; //parents ordered alphabetically
			X_ovuleNumPred[n, 2] = ovaryLen_pred[n]; //parents ordered alphabetically
			X_ovuleNumPred[n, 3] = timeToFlower_pred[n]; //parents ordered alphabetically
		}

		X_ovuleNumAll = append_array(X_ovuleNum, X_ovuleNumPred);

		f_ovuleNum_pred = gp_pred_rng(X_ovuleNumAll,ovuleNum, X_ovuleNum,
					Alpha_ovuleNum, Rho_ovuleNum, Sigma_ovuleNum, delta);

		for (n in 1:N+NPred) {
			p_ovuleNum_pred[n] = normal_rng(f_ovuleNum_pred[n], Sigma_ovuleNum);
		}

		//Split the _rep and _pred predictions
		ovuleNum_rep = p_ovuleNum_pred[1:N];
		ovuleNum_pred = p_ovuleNum_pred[N+1:N+NPred];
	}

	{
		// gynLen
		vector[NPred+N] f_gynLen_pred;
		vector[NPred+N] p_gynLen_pred;
		vector[D_gynLen] X_gynLenPred[NPred]; // declare X array for prediction
		vector[D_gynLen] X_gynLenAll[NPred+N]; // declare X array for prediction (_pred and _rep)
		// fill in X array
		for (n in 1:NPred) {
			X_gynLenPred[n, 1] = numFlowers_pred[n]; //parents ordered alphabetically
			X_gynLenPred[n, 2] = ovaryLen_pred[n]; //parents ordered alphabetically
			X_gynLenPred[n, 3] = styleLen_4pred[n]; //parents ordered alphabetically
		}

		X_gynLenAll = append_array(X_gynLen, X_gynLenPred);

		f_gynLen_pred = gp_pred_rng(X_gynLenAll,gynLen, X_gynLen,
					Alpha_gynLen, Rho_gynLen, Sigma_gynLen, delta);

		for (n in 1:N+NPred) {
			p_gynLen_pred[n] = normal_rng(f_gynLen_pred[n], Sigma_gynLen);
		}

		//Split the _rep and _pred predictions
		gynLen_rep = p_gynLen_pred[1:N];
		gynLen_pred = p_gynLen_pred[N+1:N+NPred];
	}

	{
		// percAbortedMain
		vector[NPred+N] f_percAbortedMain_pred;
		vector[NPred+N] p_percAbortedMain_pred;
		vector[D_percAbortedMain] X_percAbortedMainPred[NPred]; // declare X array for prediction
		vector[D_percAbortedMain] X_percAbortedMainAll[NPred+N]; // declare X array for prediction (_pred and _rep)
		// fill in X array
		for (n in 1:NPred) {
			X_percAbortedMainPred[n, 1] = numFlowers_pred[n]; //parents ordered alphabetically
		}

		X_percAbortedMainAll = append_array(X_percAbortedMain, X_percAbortedMainPred);

		f_percAbortedMain_pred = gp_pred_rng(X_percAbortedMainAll,percAbortedMain, X_percAbortedMain,
					Alpha_percAbortedMain, Rho_percAbortedMain, Sigma_percAbortedMain, delta);

		for (n in 1:N+NPred) {
			p_percAbortedMain_pred[n] = normal_rng(f_percAbortedMain_pred[n], Sigma_percAbortedMain);
		}

		//Split the _rep and _pred predictions
		percAbortedMain_rep = p_percAbortedMain_pred[1:N];
		percAbortedMain_pred = p_percAbortedMain_pred[N+1:N+NPred];
	}

	{
		// ovuleArea
		vector[NPred+N] f_ovuleArea_pred;
		vector[NPred+N] p_ovuleArea_pred;
		vector[D_ovuleArea] X_ovuleAreaPred[NPred]; // declare X array for prediction
		vector[D_ovuleArea] X_ovuleAreaAll[NPred+N]; // declare X array for prediction (_pred and _rep)
		// fill in X array
		for (n in 1:NPred) {
			X_ovuleAreaPred[n, 1] = gynLen_pred[n]; //parents ordered alphabetically
		}

		X_ovuleAreaAll = append_array(X_ovuleArea, X_ovuleAreaPred);

		f_ovuleArea_pred = gp_pred_rng(X_ovuleAreaAll,ovuleArea, X_ovuleArea,
					Alpha_ovuleArea, Rho_ovuleArea, Sigma_ovuleArea, delta);

		for (n in 1:N+NPred) {
			p_ovuleArea_pred[n] = normal_rng(f_ovuleArea_pred[n], Sigma_ovuleArea);
		}

		//Split the _rep and _pred predictions
		ovuleArea_rep = p_ovuleArea_pred[1:N];
		ovuleArea_pred = p_ovuleArea_pred[N+1:N+NPred];
	}

	{
		// numPodsMain
		vector[NPred+N] f_numPodsMain_pred;
		vector[NPred+N] p_numPodsMain_pred;
		vector[D_numPodsMain] X_numPodsMainPred[NPred]; // declare X array for prediction
		vector[D_numPodsMain] X_numPodsMainAll[NPred+N]; // declare X array for prediction (_pred and _rep)
		// fill in X array
		for (n in 1:NPred) {
			X_numPodsMainPred[n, 1] = num2aryBranch_4pred[n]; //parents ordered alphabetically
			X_numPodsMainPred[n, 2] = numFlowers_pred[n]; //parents ordered alphabetically
			X_numPodsMainPred[n, 3] = ovuleACov_4pred[n]; //parents ordered alphabetically
			X_numPodsMainPred[n, 4] = percAbortedMain_pred[n]; //parents ordered alphabetically
		}

		X_numPodsMainAll = append_array(X_numPodsMain, X_numPodsMainPred);

		f_numPodsMain_pred = gp_pred_rng(X_numPodsMainAll,numPodsMain, X_numPodsMain,
					Alpha_numPodsMain, Rho_numPodsMain, Sigma_numPodsMain, delta);

		for (n in 1:N+NPred) {
			p_numPodsMain_pred[n] = normal_rng(f_numPodsMain_pred[n], Sigma_numPodsMain);
		}

		//Split the _rep and _pred predictions
		numPodsMain_rep = p_numPodsMain_pred[1:N];
		numPodsMain_pred = p_numPodsMain_pred[N+1:N+NPred];
	}

	{
		// podLen
		vector[NPred+N] f_podLen_pred;
		vector[NPred+N] p_podLen_pred;
		vector[D_podLen] X_podLenPred[NPred]; // declare X array for prediction
		vector[D_podLen] X_podLenAll[NPred+N]; // declare X array for prediction (_pred and _rep)
		// fill in X array
		for (n in 1:NPred) {
			X_podLenPred[n, 1] = numFlowers_pred[n]; //parents ordered alphabetically
			X_podLenPred[n, 2] = ovuleACov_4pred[n]; //parents ordered alphabetically
			X_podLenPred[n, 3] = percAbortedMain_pred[n]; //parents ordered alphabetically
			X_podLenPred[n, 4] = styleLen_4pred[n]; //parents ordered alphabetically
		}

		X_podLenAll = append_array(X_podLen, X_podLenPred);

		f_podLen_pred = gp_pred_rng(X_podLenAll,podLen, X_podLen,
					Alpha_podLen, Rho_podLen, Sigma_podLen, delta);

		for (n in 1:N+NPred) {
			p_podLen_pred[n] = normal_rng(f_podLen_pred[n], Sigma_podLen);
		}

		//Split the _rep and _pred predictions
		podLen_rep = p_podLen_pred[1:N];
		podLen_pred = p_podLen_pred[N+1:N+NPred];
	}

	{
		// percAborted2ary
		vector[NPred+N] f_percAborted2ary_pred;
		vector[NPred+N] p_percAborted2ary_pred;
		vector[D_percAborted2ary] X_percAborted2aryPred[NPred]; // declare X array for prediction
		vector[D_percAborted2ary] X_percAborted2aryAll[NPred+N]; // declare X array for prediction (_pred and _rep)
		// fill in X array
		for (n in 1:NPred) {
			X_percAborted2aryPred[n, 1] = numFlowers_pred[n]; //parents ordered alphabetically
			X_percAborted2aryPred[n, 2] = percAbortedMain_pred[n]; //parents ordered alphabetically
		}

		X_percAborted2aryAll = append_array(X_percAborted2ary, X_percAborted2aryPred);

		f_percAborted2ary_pred = gp_pred_rng(X_percAborted2aryAll,percAborted2ary, X_percAborted2ary,
					Alpha_percAborted2ary, Rho_percAborted2ary, Sigma_percAborted2ary, delta);

		for (n in 1:N+NPred) {
			p_percAborted2ary_pred[n] = normal_rng(f_percAborted2ary_pred[n], Sigma_percAborted2ary);
		}

		//Split the _rep and _pred predictions
		percAborted2ary_rep = p_percAborted2ary_pred[1:N];
		percAborted2ary_pred = p_percAborted2ary_pred[N+1:N+NPred];
	}

	{
		// beakLen
		vector[NPred+N] f_beakLen_pred;
		vector[NPred+N] p_beakLen_pred;
		vector[D_beakLen] X_beakLenPred[NPred]; // declare X array for prediction
		vector[D_beakLen] X_beakLenAll[NPred+N]; // declare X array for prediction (_pred and _rep)
		// fill in X array
		for (n in 1:NPred) {
			X_beakLenPred[n, 1] = numPodsMain_pred[n]; //parents ordered alphabetically
			X_beakLenPred[n, 2] = podLen_pred[n]; //parents ordered alphabetically
			X_beakLenPred[n, 3] = styleLen_4pred[n]; //parents ordered alphabetically
		}

		X_beakLenAll = append_array(X_beakLen, X_beakLenPred);

		f_beakLen_pred = gp_pred_rng(X_beakLenAll,beakLen, X_beakLen,
					Alpha_beakLen, Rho_beakLen, Sigma_beakLen, delta);

		for (n in 1:N+NPred) {
			p_beakLen_pred[n] = normal_rng(f_beakLen_pred[n], Sigma_beakLen);
		}

		//Split the _rep and _pred predictions
		beakLen_rep = p_beakLen_pred[1:N];
		beakLen_pred = p_beakLen_pred[N+1:N+NPred];
	}

	{
		// numPods2ary
		vector[NPred+N] f_numPods2ary_pred;
		vector[NPred+N] p_numPods2ary_pred;
		vector[D_numPods2ary] X_numPods2aryPred[NPred]; // declare X array for prediction
		vector[D_numPods2ary] X_numPods2aryAll[NPred+N]; // declare X array for prediction (_pred and _rep)
		// fill in X array
		for (n in 1:NPred) {
			X_numPods2aryPred[n, 1] = height_4pred[n]; //parents ordered alphabetically
			X_numPods2aryPred[n, 2] = numFlowers_pred[n]; //parents ordered alphabetically
			X_numPods2aryPred[n, 3] = ovuleArea_pred[n]; //parents ordered alphabetically
			X_numPods2aryPred[n, 4] = percAborted2ary_pred[n]; //parents ordered alphabetically
			X_numPods2aryPred[n, 5] = percAbortedMain_pred[n]; //parents ordered alphabetically
		}

		X_numPods2aryAll = append_array(X_numPods2ary, X_numPods2aryPred);

		f_numPods2ary_pred = gp_pred_rng(X_numPods2aryAll,numPods2ary, X_numPods2ary,
					Alpha_numPods2ary, Rho_numPods2ary, Sigma_numPods2ary, delta);

		for (n in 1:N+NPred) {
			p_numPods2ary_pred[n] = normal_rng(f_numPods2ary_pred[n], Sigma_numPods2ary);
		}

		//Split the _rep and _pred predictions
		numPods2ary_rep = p_numPods2ary_pred[1:N];
		numPods2ary_pred = p_numPods2ary_pred[N+1:N+NPred];
	}

	{
		// timeToMature
		vector[NPred+N] f_timeToMature_pred;
		vector[NPred+N] p_timeToMature_pred;
		vector[D_timeToMature] X_timeToMaturePred[NPred]; // declare X array for prediction
		vector[D_timeToMature] X_timeToMatureAll[NPred+N]; // declare X array for prediction (_pred and _rep)
		// fill in X array
		for (n in 1:NPred) {
			X_timeToMaturePred[n, 1] = percAborted2ary_pred[n]; //parents ordered alphabetically
		}

		X_timeToMatureAll = append_array(X_timeToMature, X_timeToMaturePred);

		f_timeToMature_pred = gp_pred_rng(X_timeToMatureAll,timeToMature, X_timeToMature,
					Alpha_timeToMature, Rho_timeToMature, Sigma_timeToMature, delta);

		for (n in 1:N+NPred) {
			p_timeToMature_pred[n] = normal_rng(f_timeToMature_pred[n], Sigma_timeToMature);
		}

		//Split the _rep and _pred predictions
		timeToMature_rep = p_timeToMature_pred[1:N];
		timeToMature_pred = p_timeToMature_pred[N+1:N+NPred];
	}

	{
		// tenPodNumber
		vector[NPred+N] f_tenPodNumber_pred;
		vector[NPred+N] p_tenPodNumber_pred;
		vector[D_tenPodNumber] X_tenPodNumberPred[NPred]; // declare X array for prediction
		vector[D_tenPodNumber] X_tenPodNumberAll[NPred+N]; // declare X array for prediction (_pred and _rep)
		// fill in X array
		for (n in 1:NPred) {
			X_tenPodNumberPred[n, 1] = beakLen_pred[n]; //parents ordered alphabetically
			X_tenPodNumberPred[n, 2] = numPodsMain_pred[n]; //parents ordered alphabetically
			X_tenPodNumberPred[n, 3] = ovuleNum_pred[n]; //parents ordered alphabetically
			X_tenPodNumberPred[n, 4] = podLen_pred[n]; //parents ordered alphabetically
		}

		X_tenPodNumberAll = append_array(X_tenPodNumber, X_tenPodNumberPred);

		f_tenPodNumber_pred = gp_pred_rng(X_tenPodNumberAll,tenPodNumber, X_tenPodNumber,
					Alpha_tenPodNumber, Rho_tenPodNumber, Sigma_tenPodNumber, delta);

		for (n in 1:N+NPred) {
			p_tenPodNumber_pred[n] = normal_rng(f_tenPodNumber_pred[n], Sigma_tenPodNumber);
		}

		//Split the _rep and _pred predictions
		tenPodNumber_rep = p_tenPodNumber_pred[1:N];
		tenPodNumber_pred = p_tenPodNumber_pred[N+1:N+NPred];
	}

	{
		// tenPodArea
		vector[NPred+N] f_tenPodArea_pred;
		vector[NPred+N] p_tenPodArea_pred;
		vector[D_tenPodArea] X_tenPodAreaPred[NPred]; // declare X array for prediction
		vector[D_tenPodArea] X_tenPodAreaAll[NPred+N]; // declare X array for prediction (_pred and _rep)
		// fill in X array
		for (n in 1:NPred) {
			X_tenPodAreaPred[n, 1] = beakLen_pred[n]; //parents ordered alphabetically
			X_tenPodAreaPred[n, 2] = tenPodNumber_pred[n]; //parents ordered alphabetically
		}

		X_tenPodAreaAll = append_array(X_tenPodArea, X_tenPodAreaPred);

		f_tenPodArea_pred = gp_pred_rng(X_tenPodAreaAll,tenPodArea, X_tenPodArea,
					Alpha_tenPodArea, Rho_tenPodArea, Sigma_tenPodArea, delta);

		for (n in 1:N+NPred) {
			p_tenPodArea_pred[n] = normal_rng(f_tenPodArea_pred[n], Sigma_tenPodArea);
		}

		//Split the _rep and _pred predictions
		tenPodArea_rep = p_tenPodArea_pred[1:N];
		tenPodArea_pred = p_tenPodArea_pred[N+1:N+NPred];
	}

	{
		// seedACov
		vector[NPred+N] f_seedACov_pred;
		vector[NPred+N] p_seedACov_pred;
		vector[D_seedACov] X_seedACovPred[NPred]; // declare X array for prediction
		vector[D_seedACov] X_seedACovAll[NPred+N]; // declare X array for prediction (_pred and _rep)
		// fill in X array
		for (n in 1:NPred) {
			X_seedACovPred[n, 1] = numPodsMain_pred[n]; //parents ordered alphabetically
			X_seedACovPred[n, 2] = tenPodNumber_pred[n]; //parents ordered alphabetically
		}

		X_seedACovAll = append_array(X_seedACov, X_seedACovPred);

		f_seedACov_pred = gp_pred_rng(X_seedACovAll,seedACov, X_seedACov,
					Alpha_seedACov, Rho_seedACov, Sigma_seedACov, delta);

		for (n in 1:N+NPred) {
			p_seedACov_pred[n] = normal_rng(f_seedACov_pred[n], Sigma_seedACov);
		}

		//Split the _rep and _pred predictions
		seedACov_rep = p_seedACov_pred[1:N];
		seedACov_pred = p_seedACov_pred[N+1:N+NPred];
	}

	{
		// tenPodWeight
		vector[NPred+N] f_tenPodWeight_pred;
		vector[NPred+N] p_tenPodWeight_pred;
		vector[D_tenPodWeight] X_tenPodWeightPred[NPred]; // declare X array for prediction
		vector[D_tenPodWeight] X_tenPodWeightAll[NPred+N]; // declare X array for prediction (_pred and _rep)
		// fill in X array
		for (n in 1:NPred) {
			X_tenPodWeightPred[n, 1] = tenPodArea_pred[n]; //parents ordered alphabetically
			X_tenPodWeightPred[n, 2] = tenPodNumber_pred[n]; //parents ordered alphabetically
		}

		X_tenPodWeightAll = append_array(X_tenPodWeight, X_tenPodWeightPred);

		f_tenPodWeight_pred = gp_pred_rng(X_tenPodWeightAll,tenPodWeight, X_tenPodWeight,
					Alpha_tenPodWeight, Rho_tenPodWeight, Sigma_tenPodWeight, delta);

		for (n in 1:N+NPred) {
			p_tenPodWeight_pred[n] = normal_rng(f_tenPodWeight_pred[n], Sigma_tenPodWeight);
		}

		//Split the _rep and _pred predictions
		tenPodWeight_rep = p_tenPodWeight_pred[1:N];
		tenPodWeight_pred = p_tenPodWeight_pred[N+1:N+NPred];
	}

	{
		// totSeedArea
		vector[NPred+N] f_totSeedArea_pred;
		vector[NPred+N] p_totSeedArea_pred;
		vector[D_totSeedArea] X_totSeedAreaPred[NPred]; // declare X array for prediction
		vector[D_totSeedArea] X_totSeedAreaAll[NPred+N]; // declare X array for prediction (_pred and _rep)
		// fill in X array
		for (n in 1:NPred) {
			X_totSeedAreaPred[n, 1] = tenPodArea_pred[n]; //parents ordered alphabetically
		}

		X_totSeedAreaAll = append_array(X_totSeedArea, X_totSeedAreaPred);

		f_totSeedArea_pred = gp_pred_rng(X_totSeedAreaAll,totSeedArea, X_totSeedArea,
					Alpha_totSeedArea, Rho_totSeedArea, Sigma_totSeedArea, delta);

		for (n in 1:N+NPred) {
			p_totSeedArea_pred[n] = normal_rng(f_totSeedArea_pred[n], Sigma_totSeedArea);
		}

		//Split the _rep and _pred predictions
		totSeedArea_rep = p_totSeedArea_pred[1:N];
		totSeedArea_pred = p_totSeedArea_pred[N+1:N+NPred];
	}

	{
		// totCompactness
		vector[NPred+N] f_totCompactness_pred;
		vector[NPred+N] p_totCompactness_pred;
		vector[D_totCompactness] X_totCompactnessPred[NPred]; // declare X array for prediction
		vector[D_totCompactness] X_totCompactnessAll[NPred+N]; // declare X array for prediction (_pred and _rep)
		// fill in X array
		for (n in 1:NPred) {
			X_totCompactnessPred[n, 1] = height_4pred[n]; //parents ordered alphabetically
			X_totCompactnessPred[n, 2] = seedACov_pred[n]; //parents ordered alphabetically
		}

		X_totCompactnessAll = append_array(X_totCompactness, X_totCompactnessPred);

		f_totCompactness_pred = gp_pred_rng(X_totCompactnessAll,totCompactness, X_totCompactness,
					Alpha_totCompactness, Rho_totCompactness, Sigma_totCompactness, delta);

		for (n in 1:N+NPred) {
			p_totCompactness_pred[n] = normal_rng(f_totCompactness_pred[n], Sigma_totCompactness);
		}

		//Split the _rep and _pred predictions
		totCompactness_rep = p_totCompactness_pred[1:N];
		totCompactness_pred = p_totCompactness_pred[N+1:N+NPred];
	}

	{
		// totSeedNum
		vector[NPred+N] f_totSeedNum_pred;
		vector[NPred+N] p_totSeedNum_pred;
		vector[D_totSeedNum] X_totSeedNumPred[NPred]; // declare X array for prediction
		vector[D_totSeedNum] X_totSeedNumAll[NPred+N]; // declare X array for prediction (_pred and _rep)
		// fill in X array
		for (n in 1:NPred) {
			X_totSeedNumPred[n, 1] = num2aryBranch_4pred[n]; //parents ordered alphabetically
			X_totSeedNumPred[n, 2] = numPods2ary_pred[n]; //parents ordered alphabetically
			X_totSeedNumPred[n, 3] = seedACov_pred[n]; //parents ordered alphabetically
			X_totSeedNumPred[n, 4] = tenPodNumber_pred[n]; //parents ordered alphabetically
			X_totSeedNumPred[n, 5] = totSeedArea_pred[n]; //parents ordered alphabetically
		}

		X_totSeedNumAll = append_array(X_totSeedNum, X_totSeedNumPred);

		f_totSeedNum_pred = gp_pred_rng(X_totSeedNumAll,totSeedNum, X_totSeedNum,
					Alpha_totSeedNum, Rho_totSeedNum, Sigma_totSeedNum, delta);

		for (n in 1:N+NPred) {
			p_totSeedNum_pred[n] = normal_rng(f_totSeedNum_pred[n], Sigma_totSeedNum);
		}

		//Split the _rep and _pred predictions
		totSeedNum_rep = p_totSeedNum_pred[1:N];
		totSeedNum_pred = p_totSeedNum_pred[N+1:N+NPred];
	}

	{
		// tenPodCompactness
		vector[NPred+N] f_tenPodCompactness_pred;
		vector[NPred+N] p_tenPodCompactness_pred;
		vector[D_tenPodCompactness] X_tenPodCompactnessPred[NPred]; // declare X array for prediction
		vector[D_tenPodCompactness] X_tenPodCompactnessAll[NPred+N]; // declare X array for prediction (_pred and _rep)
		// fill in X array
		for (n in 1:NPred) {
			X_tenPodCompactnessPred[n, 1] = seedACov_pred[n]; //parents ordered alphabetically
			X_tenPodCompactnessPred[n, 2] = totCompactness_pred[n]; //parents ordered alphabetically
		}

		X_tenPodCompactnessAll = append_array(X_tenPodCompactness, X_tenPodCompactnessPred);

		f_tenPodCompactness_pred = gp_pred_rng(X_tenPodCompactnessAll,tenPodCompactness, X_tenPodCompactness,
					Alpha_tenPodCompactness, Rho_tenPodCompactness, Sigma_tenPodCompactness, delta);

		for (n in 1:N+NPred) {
			p_tenPodCompactness_pred[n] = normal_rng(f_tenPodCompactness_pred[n], Sigma_tenPodCompactness);
		}

		//Split the _rep and _pred predictions
		tenPodCompactness_rep = p_tenPodCompactness_pred[1:N];
		tenPodCompactness_pred = p_tenPodCompactness_pred[N+1:N+NPred];
	}

	{
		// TGW
		vector[NPred+N] f_TGW_pred;
		vector[NPred+N] p_TGW_pred;
		vector[D_TGW] X_TGWPred[NPred]; // declare X array for prediction
		vector[D_TGW] X_TGWAll[NPred+N]; // declare X array for prediction (_pred and _rep)
		// fill in X array
		for (n in 1:NPred) {
			X_TGWPred[n, 1] = seedACov_pred[n]; //parents ordered alphabetically
			X_TGWPred[n, 2] = styleLen_4pred[n]; //parents ordered alphabetically
			X_TGWPred[n, 3] = tenPodArea_pred[n]; //parents ordered alphabetically
			X_TGWPred[n, 4] = totSeedNum_pred[n]; //parents ordered alphabetically
		}

		X_TGWAll = append_array(X_TGW, X_TGWPred);

		f_TGW_pred = gp_pred_rng(X_TGWAll,TGW, X_TGW,
					Alpha_TGW, Rho_TGW, Sigma_TGW, delta);

		for (n in 1:N+NPred) {
			p_TGW_pred[n] = normal_rng(f_TGW_pred[n], Sigma_TGW);
		}

		//Split the _rep and _pred predictions
		TGW_rep = p_TGW_pred[1:N];
		TGW_pred = p_TGW_pred[N+1:N+NPred];
	}

	{
		// totSeedW
		vector[NPred+N] f_totSeedW_pred;
		vector[NPred+N] p_totSeedW_pred;
		vector[D_totSeedW] X_totSeedWPred[NPred]; // declare X array for prediction
		vector[D_totSeedW] X_totSeedWAll[NPred+N]; // declare X array for prediction (_pred and _rep)
		// fill in X array
		for (n in 1:NPred) {
			X_totSeedWPred[n, 1] = numFlowers_pred[n]; //parents ordered alphabetically
			X_totSeedWPred[n, 2] = numPods2ary_pred[n]; //parents ordered alphabetically
			X_totSeedWPred[n, 3] = numPodsMain_pred[n]; //parents ordered alphabetically
			X_totSeedWPred[n, 4] = tenPodWeight_pred[n]; //parents ordered alphabetically
			X_totSeedWPred[n, 5] = TGW_pred[n]; //parents ordered alphabetically
			X_totSeedWPred[n, 6] = totSeedArea_pred[n]; //parents ordered alphabetically
			X_totSeedWPred[n, 7] = totSeedNum_pred[n]; //parents ordered alphabetically
		}

		X_totSeedWAll = append_array(X_totSeedW, X_totSeedWPred);

		f_totSeedW_pred = gp_pred_rng(X_totSeedWAll,totSeedW, X_totSeedW,
					Alpha_totSeedW, Rho_totSeedW, Sigma_totSeedW, delta);

		for (n in 1:N+NPred) {
			p_totSeedW_pred[n] = normal_rng(f_totSeedW_pred[n], Sigma_totSeedW);
		}

		//Split the _rep and _pred predictions
		totSeedW_rep = p_totSeedW_pred[1:N];
		totSeedW_pred = p_totSeedW_pred[N+1:N+NPred];
	}

	{
		// oilContent
		vector[NPred+N] f_oilContent_pred;
		vector[NPred+N] p_oilContent_pred;
		vector[D_oilContent] X_oilContentPred[NPred]; // declare X array for prediction
		vector[D_oilContent] X_oilContentAll[NPred+N]; // declare X array for prediction (_pred and _rep)
		// fill in X array
		for (n in 1:NPred) {
			X_oilContentPred[n, 1] = ovuleNum_pred[n]; //parents ordered alphabetically
			X_oilContentPred[n, 2] = totCompactness_pred[n]; //parents ordered alphabetically
			X_oilContentPred[n, 3] = totSeedW_pred[n]; //parents ordered alphabetically
		}

		X_oilContentAll = append_array(X_oilContent, X_oilContentPred);

		f_oilContent_pred = gp_pred_rng(X_oilContentAll,oilContent, X_oilContent,
					Alpha_oilContent, Rho_oilContent, Sigma_oilContent, delta);

		for (n in 1:N+NPred) {
			p_oilContent_pred[n] = normal_rng(f_oilContent_pred[n], Sigma_oilContent);
		}

		//Split the _rep and _pred predictions
		oilContent_rep = p_oilContent_pred[1:N];
		oilContent_pred = p_oilContent_pred[N+1:N+NPred];
	}

}

